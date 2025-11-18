// Copyright 2019-2025 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <TRandom.h>
#include <cmath>
#include <numeric>
#include "FD3Simulation/Digitizer.h"
#include "FD3Simulation/DigitizationConstants.h"
#include "FD3Simulation/FD3DigParam.h"
#include "FD3Base/GeometryTGeo.h"
#include "FD3Base/Constants.h"
#include "TF1Convolution.h"

ClassImp(o2::fd3::Digitizer);

using namespace o2::math_utils;
using namespace o2::fd3;

void Digitizer::clear()
{
  mEventId = -1;
  mSrcId = -1;
  for (auto& analogSignal : mPmtChargeVsTime) {
    std::fill_n(std::begin(analogSignal), analogSignal.size(), 0);
  }
  mLastBCCache.clear();
  mCfdStartIndex.fill(0);
}

//_______________________________________________________________________
void Digitizer::init()
{
  LOG(info) << "init";
  mNBins = FD3DigParam::Instance().waveformNbins;      //Will be computed using detector set-up from CDB
  mBinSize = FD3DigParam::Instance().waveformBinWidth; //Will be set-up from CDB
  mNTimeBinsPerBC = std::lround(o2::constants::lhc::LHCBunchSpacingNS / mBinSize); // 1920 bins/BC

  for (Int_t detID = 0; detID < DigitizationConstants::NCELLSTOT; detID++) {
    mPmtChargeVsTime[detID].resize(mNBins);
    mLastBCCache.mPmtChargeVsTime[detID].resize(mNBins);
  }

  /// set up PMT response function [avg] for all rings
  TF1Convolution convolutionRings("expo", "landau", 5.e-09, 90.e-09, false);
  TF1 convolutionRingsFn("convolutionFn", convolutionRings, 5.e-09, 90.e-09, convolutionRings.GetNpar());
  convolutionRingsFn.SetParameters(FD3DigParam::Instance().constRings, FD3DigParam::Instance().slopeRings,
                                        FD3DigParam::Instance().mpvRings, FD3DigParam::Instance().sigmaRings);

  /// PMT response per hit [Global] for ring 1 to 4
  mPmtResponseGlobalRings.resize(mNBins);
  const float binSizeInNs = mBinSize * 1.e-09; // to convert ns into sec
  double x = (binSizeInNs) / 2.0;
  for (auto& y : mPmtResponseGlobalRings) {
    y = FD3DigParam::Instance().getNormRings()                                   // normalisation to have MIP adc at 16
        * convolutionRingsFn.Eval(x + FD3DigParam::Instance().offsetRings); // offset to adjust mean position of waveform
    x += binSizeInNs;
  }

  mLastBCCache.clear();
  mCfdStartIndex.fill(0);

  LOG(info) << "init -> finished";
}

//_______________________________________________________________________
void Digitizer::process(const std::vector<o2::fd3::Hit>& hits,
                        std::vector<o2::fd3::Digit>& digitsBC,
                        std::vector<o2::fd3::ChannelData>& digitsCh,
                        std::vector<o2::fd3::DetTrigInput>& digitsTrig,
                        o2::dataformats::MCTruthContainer<o2::fd3::MCLabel>& labels)
{
  LOG(debug) << "Begin with " << hits.size() << " hits";
  flush(digitsBC, digitsCh, digitsTrig, labels); // flush cached signal which cannot be affect by new event

  std::vector<int> hitIdx(hits.size());
  std::iota(std::begin(hitIdx), std::end(hitIdx), 0);
  std::sort(std::begin(hitIdx), std::end(hitIdx),
            [&hits](int a, int b) { return hits[a].GetTrackID() < hits[b].GetTrackID(); });

  // use ordered hits
  for (auto ids : hitIdx) {
    const auto& hit = hits[ids];
    Int_t detId = hit.GetDetectorID();
    double hitEdep = hit.GetHitValue() * 1e3;  // convert to MeV
    float const hitTime = hit.GetTime() * 1e9; // convert to ns
    // TODO: check how big is inaccuracy if more than 1 'below-threshold' particles hit the same detector cell
    if (hitEdep < FD3DigParam::Instance().singleMipThreshold || hitTime > FD3DigParam::Instance().singleHitTimeThreshold) {
      continue;
    }
    float distanceFromXc = 0;

    double const nPhotons = hitEdep * DP::N_PHOTONS_PER_MEV;
    float const nPhE = SimulateLightYield(detId, nPhotons);
    float const mipFraction = float(nPhE / FD3DigParam::Instance().avgNumberPhElectronPerMip);
    Long64_t timeHit = hitTime;
    timeHit += mIntRecord.getTimeNS();
    o2::InteractionTimeRecord const irHit(timeHit);
    std::array<o2::InteractionRecord, NBC2Cache> cachedIR;
    int nCachedIR = 0;
    for (int i = BCCacheMin; i < BCCacheMax + 1; i++) {
      double const tNS = timeHit + o2::constants::lhc::LHCBunchSpacingNS * i;
      cachedIR[nCachedIR].setFromNS(tNS);
      if (tNS < 0 && cachedIR[nCachedIR] > irHit) {
        continue; // don't go to negative BC/orbit (it will wrap)
      }
      // ensure existence of cached container
      setBCCache(cachedIR[nCachedIR++]);
    } // BCCache loop

    createPulse(mipFraction, hit.GetTrackID(), hitTime, hit.GetPos().R(), cachedIR, nCachedIR, detId);

  }   //hitloop
}

//_______________________________________________________________________
void Digitizer::createPulse(float mipFraction, int parID, const double hitTime, const float hitR,
                            std::array<o2::InteractionRecord, NBC2Cache> const& cachedIR, int nCachedIR, const int detId)
{

  std::array<bool, NBC2Cache> added;
  added.fill(false);

  for (int ir = 0; ir < NBC2Cache; ir++) {
    auto bcCache = getBCCache(cachedIR[ir]);
    for (int ich = 0; ich < DigitizationConstants::NCELLSTOT; ich++) {
      (*bcCache).mPmtChargeVsTime[ich].resize(mNTimeBinsPerBC);
    }
  }

  // Subtract time-of-flight from hit time
  const float timeOfFlight = hitR / o2::constants::physics::LightSpeedCm2NS;
  Int_t const NBinShift = std::lround((hitTime - timeOfFlight + FD3DigParam::Instance().hitTimeOffset) / FD3DigParam::Instance().waveformBinWidth);

  if (NBinShift >= 0 && NBinShift < FD3DigParam::Instance().waveformNbins) {
    mPmtResponseTemp.resize(FD3DigParam::Instance().waveformNbins, 0.);
    std::memcpy(&mPmtResponseTemp[NBinShift], &mPmtResponseGlobalRings[0],
              sizeof(double) * (FD3DigParam::Instance().waveformNbins - NBinShift));
  } else {
    mPmtResponseTemp = mPmtResponseGlobalRings;
    mPmtResponseTemp.erase(mPmtResponseTemp.begin(), mPmtResponseTemp.begin() + abs(NBinShift));
    mPmtResponseTemp.resize(FD3DigParam::Instance().waveformNbins);
  }

  for (int ir = 0; ir < int(mPmtResponseTemp.size() / mNTimeBinsPerBC); ir++) {
    auto bcCache = getBCCache(cachedIR[ir]);

    for (int iBin = 0; iBin < mNTimeBinsPerBC; iBin++) {
      (*bcCache).mPmtChargeVsTime[detId][iBin] += (mPmtResponseTemp[ir * mNTimeBinsPerBC + iBin] * mipFraction);
    }
    added[ir] = true;
  }
  // Add MC labels to BCs for those contributed to the PMT signal
  for (int ir = 0; ir < nCachedIR; ir++) {
    if (added[ir]) {
      auto bcCache = getBCCache(cachedIR[ir]);
      (*bcCache).labels.emplace_back(parID, mEventId, mSrcId, detId);
    }
  }
}

//_______________________________________________________________________
void Digitizer::flush(std::vector<o2::fd3::Digit>& digitsBC,
                      std::vector<o2::fd3::ChannelData>& digitsCh,
                      std::vector<o2::fd3::DetTrigInput>& digitsTrig,
                      o2::dataformats::MCTruthContainer<o2::fd3::MCLabel>& labels)
{
  ++mEventId;
  while (!mCache.empty()) {
    auto const& bc = mCache.front();
    if (mIntRecord.differenceInBC(bc) > NBC2Cache) { // Build events that are separated by NBC2Cache BCs from current BC
      storeBC(bc, digitsBC, digitsCh, digitsTrig, labels);
      mCache.pop_front();
    } else {
      return;
    }
  }
}

//_______________________________________________________________________
void Digitizer::storeBC(const BCCache& bc,
                        std::vector<o2::fd3::Digit>& digitsBC,
                        std::vector<o2::fd3::ChannelData>& digitsCh,
                        std::vector<o2::fd3::DetTrigInput>& digitsTrig,
                        o2::dataformats::MCTruthContainer<o2::fd3::MCLabel>& labels)

{
  size_t const nBC = digitsBC.size();   // save before digitsBC is being modified
  size_t const first = digitsCh.size(); // save before digitsCh is being modified
  int8_t nTotFiredCells = 0;
  int8_t nTrgFiredCells = 0; // number of fired cells, that follow additional trigger conditions (time gate)
  int totalChargeAllRing = 0;
  int32_t avgTime = 0;
  double nSignalInner = 0;
  double nSignalOuter = 0;

  int nhitsA = 0, nhitsC = 0, meanTimeA = 0, meanTimeC = 0;
  int sumAmplA = 0, sumAmplC = 0;
  double nSignalA = 0;
  double nSignalC = 0;

  if (mLastBCCache.differenceInBC(bc) != 1) { // if the last buffered BC is not the one before the current BC
    mLastBCCache.clear();                     // clear the bufffer (mPmtChargeVsTime set to 0s)
    mCfdStartIndex.fill(0);                   // reset all start indices to 0, i.e., to the beginning of the BC
  }

  for (int iPmt = 0; iPmt < DigitizationConstants::NCELLSTOT; iPmt++) {
    // run the CFD: this updates the start index for the next BC in case the CFD dead time ends in the next BC
    double cfdWithOffset = SimulateTimeCfd(mCfdStartIndex[iPmt], mLastBCCache.mPmtChargeVsTime[iPmt], bc.mPmtChargeVsTime[iPmt]);
    double cfdZero = cfdWithOffset - FD3DigParam::Instance().avgCfdTimeForMip;

    // Conditions to sum charge are: all participating channels must have time within +/- 2.5 ns, AND
    //   at least one channel must follow more strict conditions (see below)
    if (cfdZero < -FD3DigParam::Instance().cfdCheckWindow || cfdZero > FD3DigParam::Instance().cfdCheckWindow) {
      continue;
    }

    int iTotalCharge = std::lround(IntegrateCharge(bc.mPmtChargeVsTime[iPmt]) * DP::INV_CHARGE_PER_ADC); // convert Coulomb to adc;

    uint8_t channelBits = FD3DigParam::Instance().defaultChainQtc;
    if (std::rand() % 2) {
      ChannelData::setFlag(ChannelData::kNumberADC, channelBits);
    }
    if (iTotalCharge > (FD3DigParam::Instance().maxCountInAdc) && FD3DigParam::Instance().useMaxChInAdc) {
      iTotalCharge = FD3DigParam::Instance().maxCountInAdc; // max adc channel for one PMT
      ChannelData::setFlag(ChannelData::kIsAmpHigh, channelBits);
    }

    if (iTotalCharge < FD3DigParam::Instance().getCFDTrshInAdc()) {
      continue;
    }

    int iCfdZero = std::lround(cfdZero * DP::INV_TIME_PER_TDCCHANNEL);
    digitsCh.emplace_back(iPmt, iCfdZero, iTotalCharge, channelBits);
    ++nTotFiredCells;

    int triggerGate = FD3DigParam::Instance().mTime_trg_gate;
    if (std::abs(iCfdZero) < triggerGate) {
      ++nTrgFiredCells;
      //---trigger---
      totalChargeAllRing += iTotalCharge;
      avgTime += iCfdZero;
      if (iPmt < DigitizationConstants::NCELLSA) {
         nSignalA++;
      } else {
         nSignalC++;
      }
    }
  }
  // save BC information for the CFD detector
  mLastBCCache = bc;
  if (nTotFiredCells < 1) {
    return;
  }
  if (nTrgFiredCells > 0) {
    avgTime /= nTrgFiredCells;
  } else {
    avgTime = o2::fd3::Triggers::DEFAULT_TIME;
  }
  ///Triggers for FD3
  bool isA, isAIn, isAOut, isCen, isSCen;
  isA = nTrgFiredCells > 0;
  isAIn = nSignalInner > 0;  // ring 1,2 and 3
  isAOut = nSignalOuter > 0; // ring 4 and 5
  isCen = totalChargeAllRing > FD3DigParam::Instance().adcChargeCenThr;
  isSCen = totalChargeAllRing > FD3DigParam::Instance().adcChargeSCenThr;

  Triggers triggers;
  const int unusedCharge = o2::fd3::Triggers::DEFAULT_AMP;
  const int unusedTime = o2::fd3::Triggers::DEFAULT_TIME;
  const int unusedZero = o2::fd3::Triggers::DEFAULT_ZERO;
  const bool unusedBitsInSim = false; // bits related to laser and data validity
  const bool bitDataIsValid = true;
  triggers.setTriggers(isA, isAIn, isAOut, isCen, isSCen, nTrgFiredCells, (int8_t)unusedZero,
                       (int32_t)(0.125 * totalChargeAllRing), (int32_t)unusedCharge, (int16_t)avgTime, (int16_t)unusedTime, unusedBitsInSim, unusedBitsInSim, bitDataIsValid);
  digitsBC.emplace_back(first, nTotFiredCells, bc, triggers, mEventId - 1);
  digitsTrig.emplace_back(bc, isA, isAIn, isAOut, isCen, isSCen);
  for (auto const& lbl : bc.labels) {
    labels.addElement(nBC, lbl);
  }
}

// -------------------------------------------------------------------------------
// --- Internal helper methods related to conversion of energy-deposition into ---
// --- photons -> photoelectrons -> electrical signal                          ---
// -------------------------------------------------------------------------------
Int_t Digitizer::SimulateLightYield(Int_t pmt, Int_t nPhot) const
{
  const float epsilon = 0.0001f;
  const float p = FD3DigParam::Instance().lightYield * FD3DigParam::Instance().photoCathodeEfficiency;
  if ((fabs(1.0f - p) < epsilon) || nPhot == 0) {
    return nPhot;
  }
  const Int_t n = Int_t(nPhot < 100
                          ? gRandom->Binomial(nPhot, p)
                          : gRandom->Gaus((p * nPhot) + 0.5, TMath::Sqrt(p * (1. - p) * nPhot)));
  return n;
}
//---------------------------------------------------------------------------
float Digitizer::IntegrateCharge(const ChannelDigitF& pulse) const
{
  int const chargeIntMin = FD3DigParam::Instance().isIntegrateFull ? 0 : (FD3DigParam::Instance().avgCfdTimeForMip - 6.0) / mBinSize;                //Charge integration offset (cfd mean time - 6 ns)
  int const chargeIntMax = FD3DigParam::Instance().isIntegrateFull ? mNTimeBinsPerBC : (FD3DigParam::Instance().avgCfdTimeForMip + 14.0) / mBinSize; //Charge integration offset (cfd mean time + 14 ns)
  if (chargeIntMin < 0 || chargeIntMin > mNTimeBinsPerBC || chargeIntMax > mNTimeBinsPerBC) {
    LOG(fatal) << "invalid indicess: chargeInMin=" << chargeIntMin << " chargeIntMax=" << chargeIntMax;
  }
  float totalCharge = 0.0f;
  for (int iTimeBin = chargeIntMin; iTimeBin < chargeIntMax; iTimeBin++) {
    totalCharge += pulse[iTimeBin];
  }
  return totalCharge;
}

//---------------------------------------------------------------------------
float Digitizer::SimulateTimeCfd(int& startIndex, const ChannelDigitF& pulseLast, const ChannelDigitF& pulse) const
{
  float timeCfd = -1024.0f;

  if (pulse.empty()) {
    startIndex = 0;
    return timeCfd;
  }

  float const cfdThrInCoulomb = FD3DigParam::Instance().mCFD_trsh * 1e-3 / 50 * mBinSize * 1e-9; // convert mV into Coulomb assuming 50 Ohm

  Int_t const binShift = TMath::Nint(FD3DigParam::Instance().timeShiftCfd / mBinSize);
  float sigPrev = 5 * pulseLast[mNTimeBinsPerBC - binShift - 1] - pulseLast[mNTimeBinsPerBC - 1]; //  CFD output from the last bin of the last BC
  for (Int_t iTimeBin = 0; iTimeBin < mNTimeBinsPerBC; ++iTimeBin) {
    float const sigCurrent = 5.0f * (iTimeBin >= binShift ? pulse[iTimeBin - binShift] : pulseLast[mNTimeBinsPerBC - binShift + iTimeBin]) - pulse[iTimeBin];
    if (iTimeBin >= startIndex && std::abs(pulse[iTimeBin]) > cfdThrInCoulomb) { // enable
      if (sigPrev < 0.0f && sigCurrent >= 0.0f) {                                // test for zero-crossing
        timeCfd = float(iTimeBin) * mBinSize;
        startIndex = iTimeBin + std::lround(FD3DigParam::Instance().mCfdDeadTime / mBinSize); // update startIndex (CFD dead time)
        if (startIndex < mNTimeBinsPerBC) {
          startIndex = 0; // dead-time ends in same BC: no impact on the following BC
        } else {
          startIndex -= mNTimeBinsPerBC;
        }
        if (startIndex > mNTimeBinsPerBC) {
          LOG(fatal) << "CFD dead-time was set to > 25 ns";
        }
        break; // only detects the 1st zero-crossing in the BC
      }
    }
    sigPrev = sigCurrent;
  }
  return timeCfd;
}

//---------------------------------------------------------------------------
//float Digitizer::getDistFromCellCenter(unsigned int cellId, double hitx, double hity)
//{
//  Geometry* geo = Geometry::instance();
//
//  // Parametrize the line (ax+by+c=0) that crosses the detector center and the cell's middle point
//  Point3Dsimple* pCell = &geo->getCellCenter(cellId);
//  float x0, y0, z0;
//  geo->getGlobalPosition(x0, y0, z0);
//  double a = -(y0 - pCell->y) / (x0 - pCell->x);
//  double b = 1;
//  double c = -(y0 - a * x0);
//  //Return the distance from hit to this line
//  return (a * hitx + b * hity + c) / TMath::Sqrt(a * a + b * b);
//  return 0.;
//}

//float Digitizer::getSignalFraction(float distanceFromXc, bool isFirstChannel)
//{
 // return 1.;
//  float const fraction = sigmoidPmtRing5(distanceFromXc);
//  if (distanceFromXc > 0) {
//    return isFirstChannel ? fraction : (1. - fraction);
//  } else {
//    return isFirstChannel ? (1. - fraction) : fraction;
//  }
//}

//_____________________________________________________________________________
o2::fd3::Digitizer::BCCache& Digitizer::setBCCache(const o2::InteractionRecord& ir)
{
  if (mCache.empty() || mCache.back() < ir) {
    mCache.emplace_back();
    auto& cb = mCache.back();
    cb = ir;
    return cb;
  }
  if (mCache.front() > ir) {
    mCache.emplace_front();
    auto& cb = mCache.front();
    cb = ir;
    return cb;
  }
  for (auto cb = mCache.begin(); cb != mCache.end(); cb++) {
    if ((*cb) == ir) {
      return *cb;
    }
    if (ir < (*cb)) {
      auto cbnew = mCache.emplace(cb); // insert new element before cb
      (*cbnew) = ir;
      return (*cbnew);
    }
  }
  return mCache.front();
}
//_____________________________________________________________________________
o2::fd3::Digitizer::BCCache* Digitizer::getBCCache(const o2::InteractionRecord& ir)
{
  // get pointer on existing cache
  for (auto cb = mCache.begin(); cb != mCache.end(); cb++) {
    if ((*cb) == ir) {
      return &(*cb);
    }
  }
  return nullptr;
}

//bool Digitizer::isRing5(int detID)
//{
//  if (detID > 31) {
//    return true;
//  } else {
//    return false;
//  }
//}

O2ParamImpl(FD3DigParam);
