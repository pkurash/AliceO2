// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
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
#include "FD2Simulation/DigitizationConstant.h"
#include "FD2Simulation/Digitizer.h"
#include "FD2Simulation/FD2DigParam.h"
#include "FD2Base/GeometryTGeo.h"
#include "FD2Base/Constants.h"
#include "TF1Convolution.h"

ClassImp(o2::fd2::Digitizer);

using namespace o2::math_utils;
using namespace o2::fd2;

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
  mNBins = FD2DigParam::Instance().waveformNbins;                                  // Will be computed using detector set-up from CDB
  mBinSize = FD2DigParam::Instance().waveformBinWidth;                             // Will be set-up from CDB
  mNTimeBinsPerBC = std::lround(o2::constants::lhc::LHCBunchSpacingNS / mBinSize); // 1920 bins/BC

//  for (Int_t detID = 0; detID < Constants::nFv0Channels; detID++) {
  for (Int_t detID = 0; detID <  16; detID++) {
    mPmtChargeVsTime[detID].resize(mNBins);
    mLastBCCache.mPmtChargeVsTime[detID].resize(mNBins);
  }

  /// set up PMT response function [avg] for ring 1 to 4
  TF1Convolution convolutionRingA1ToA4("expo", "landau", 5.e-09, 90.e-09, false);
  TF1 convolutionRingA1ToA4Fn("convolutionFn", convolutionRingA1ToA4, 5.e-09, 90.e-09, convolutionRingA1ToA4.GetNpar());
  convolutionRingA1ToA4Fn.SetParameters(FD2DigParam::Instance().constRingA1ToA4, FD2DigParam::Instance().slopeRingA1ToA4,
                                        FD2DigParam::Instance().mpvRingA1ToA4, FD2DigParam::Instance().sigmaRingA1ToA4);

  /// set up PMT response function [avg] for ring 5
  TF1Convolution convolutionRing5("expo", "landau", 5.e-09, 90.e-09, false);
  TF1 convolutionRing5Fn("convolutionFn", convolutionRing5, 5.e-09, 90.e-09, convolutionRing5.GetNpar());
  convolutionRing5Fn.SetParameters(FD2DigParam::Instance().constRing5, FD2DigParam::Instance().slopeRing5,
                                   FD2DigParam::Instance().mpvRing5, FD2DigParam::Instance().sigmaRing5);
  /// PMT response per hit [Global] for ring 1 to 4
  mPmtResponseGlobalRingA1ToA4.resize(mNBins);
  const float binSizeInNs = mBinSize * 1.e-09; // to convert ns into sec
  double x = (binSizeInNs) / 2.0;
  for (auto& y : mPmtResponseGlobalRingA1ToA4) {
    y = FD2DigParam::Instance().getNormRingA1ToA4()                                   // normalisation to have MIP adc at 16
        * convolutionRingA1ToA4Fn.Eval(x + FD2DigParam::Instance().offsetRingA1ToA4); // offset to adjust mean position of waveform
    x += binSizeInNs;
  }
  /// PMT response per hit [Global] for ring 5
  mPmtResponseGlobalRing5.resize(mNBins);
  x = (binSizeInNs) / 2.0;
  for (auto& y : mPmtResponseGlobalRing5) {
    y = FD2DigParam::Instance().getNormRing5()                              // normalisation to have MIP adc at 16
        * convolutionRing5Fn.Eval(x + FD2DigParam::Instance().offsetRing5); // offset to adjust mean position of waveform
    x += binSizeInNs;
  }
  mLastBCCache.clear();
  mCfdStartIndex.fill(0);

  InitCellCenters();

  LOG(info) << "init -> finished";
}

void Digitizer::process(const std::vector<o2::fd3::Hit>& hits,
                        std::vector<o2::fv0::Digit>& digitsBC,
                        std::vector<o2::fv0::ChannelData>& digitsCh,
                        std::vector<o2::fv0::DetTrigInput>& digitsTrig,
                        o2::dataformats::MCTruthContainer<o2::fv0::MCLabel>& labels)
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

    if (mDeadChannelMap && !mDeadChannelMap->isChannelAlive(detId)) {
      continue;
    }

    Double_t hitEdep = hit.GetHitValue() * 1e3;  // convert to MeV
    Float_t const hitTime = hit.GetTime() * 1e9; // convert to ns
    // TODO: check how big is inaccuracy if more than 1 'below-threshold' particles hit the same detector cell
    if (hitEdep < FD2DigParam::Instance().singleMipThreshold || hitTime > FD2DigParam::Instance().singleHitTimeThreshold) {
      continue;
    }
    float distanceFromXc = 0;
    //if (Geometry::instance()->isRing5(detId)) {
    //  distanceFromXc = getDistFromCellCenter(detId, hit.GetX(), hit.GetY());
    //}

    int iChannelPerCell = 0;
    while (iChannelPerCell < 2) { // loop over 2 channels, into which signal from each cell in ring 5 is split
    //  if (Geometry::instance()->isRing5(detId)) {
    //    // The first channel number is located counter-clockwise from the cell center
    //    // and remains identical to the detector number, the second one is clockwise and incremented by 8
    //    if (iChannelPerCell == 1) {
    //      detId += 8;
    //    }
    //    // Split signal magnitude to fractions depending on the distance of the hit from the cell center
    //    hitEdep = (hit.GetHitValue() * 1e3) * getSignalFraction(distanceFromXc, iChannelPerCell == 0);
    //    // LOG(info) << "  detId: " << detId << "-" << iChannelPerCell << " hitEdep: " << hitEdep << " distanceFromXc: " << distanceFromXc;
    //    ++iChannelPerCell;
    //  } else {
    //    iChannelPerCell = 2; // not a ring 5 cell -> don't repeat the loop
    //  }
      Double_t const nPhotons = hitEdep * DP::N_PHOTONS_PER_MEV;
      float const nPhE = SimulateLightYield(detId, nPhotons);
      float const mipFraction = float(nPhE / FD2DigParam::Instance().avgNumberPhElectronPerMip);
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

    } // while loop
  } // hitloop
}

void Digitizer::createPulse(float mipFraction, int parID, const double hitTime, const float hitR,
                            std::array<o2::InteractionRecord, NBC2Cache> const& cachedIR, int nCachedIR, const int detId)
{

  std::array<bool, NBC2Cache> added;
  added.fill(false);

  for (int ir = 0; ir < NBC2Cache; ir++) {
    auto bcCache = getBCCache(cachedIR[ir]);
//    for (int ich = 0; ich < Constants::nFv0Channels; ich++) {
    for (int ich = 0; ich < 16; ich++) {
      (*bcCache).mPmtChargeVsTime[ich].resize(mNTimeBinsPerBC);
    }
  }

  // Subtract time-of-flight from hit time
  const float timeOfFlight = hitR / o2::constants::physics::LightSpeedCm2NS;
  Int_t const NBinShift = std::lround((hitTime - timeOfFlight + FD2DigParam::Instance().hitTimeOffset) / FD2DigParam::Instance().waveformBinWidth);

  if (NBinShift >= 0 && NBinShift < FD2DigParam::Instance().waveformNbins) {
    mPmtResponseTemp.resize(FD2DigParam::Instance().waveformNbins, 0.);
    if (isRing5(detId)) {
      std::memcpy(&mPmtResponseTemp[NBinShift], &mPmtResponseGlobalRing5[0],
                  sizeof(double) * (FD2DigParam::Instance().waveformNbins - NBinShift));
    } else {
      std::memcpy(&mPmtResponseTemp[NBinShift], &mPmtResponseGlobalRingA1ToA4[0],
                  sizeof(double) * (FD2DigParam::Instance().waveformNbins - NBinShift));
    }
  } else {
    if (isRing5(detId)) {
      mPmtResponseTemp = mPmtResponseGlobalRing5;
      mPmtResponseTemp.erase(mPmtResponseTemp.begin(), mPmtResponseTemp.begin() + abs(NBinShift));
    } else {
      mPmtResponseTemp = mPmtResponseGlobalRingA1ToA4;
      mPmtResponseTemp.erase(mPmtResponseTemp.begin(), mPmtResponseTemp.begin() + abs(NBinShift));
    }

    mPmtResponseTemp.resize(FD2DigParam::Instance().waveformNbins);
  }

  for (int ir = 0; ir < int(mPmtResponseTemp.size() / mNTimeBinsPerBC); ir++) {
    auto bcCache = getBCCache(cachedIR[ir]);

    for (int iBin = 0; iBin < mNTimeBinsPerBC; iBin++) {
      (*bcCache).mPmtChargeVsTime[detId][iBin] += (mPmtResponseTemp[ir * mNTimeBinsPerBC + iBin] * mipFraction);
    }
    added[ir] = true;
  }
  /// Add MC labels to BCs for those contributed to the PMT signal
  for (int ir = 0; ir < nCachedIR; ir++) {
    if (added[ir]) {
      auto bcCache = getBCCache(cachedIR[ir]);
      (*bcCache).labels.emplace_back(parID, mEventId, mSrcId, detId);
    }
  }
}

void Digitizer::flush(std::vector<o2::fv0::Digit>& digitsBC,
                      std::vector<o2::fv0::ChannelData>& digitsCh,
                      std::vector<o2::fv0::DetTrigInput>& digitsTrig,
                      o2::dataformats::MCTruthContainer<o2::fv0::MCLabel>& labels)
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

void Digitizer::storeBC(const BCCache& bc,
                        std::vector<o2::fv0::Digit>& digitsBC,
                        std::vector<o2::fv0::ChannelData>& digitsCh,
                        std::vector<o2::fv0::DetTrigInput>& digitsTrig,
                        o2::dataformats::MCTruthContainer<o2::fv0::MCLabel>& labels)

{
  size_t const nBC = digitsBC.size();   // save before digitsBC is being modified
  size_t const first = digitsCh.size(); // save before digitsCh is being modified
  int8_t nTotFiredCells = 0;
  int8_t nTrgFiredCells = 0; // number of fired cells, that follow additional trigger conditions (time gate)
  int totalChargeAllRing = 0;
  int totalChargeInnerRing = 0;
  int totalChargeOuterRing = 0;
  int32_t avgTime = 0;
  double nSignalInner = 0;
  double nSignalOuter = 0;

  if (mLastBCCache.differenceInBC(bc) != 1) { // if the last buffered BC is not the one before the current BC
    mLastBCCache.clear();                     // clear the bufffer (mPmtChargeVsTime set to 0s)
    mCfdStartIndex.fill(0);                   // reset all start indices to 0, i.e., to the beginning of the BC
  }

//  for (int iPmt = 0; iPmt < Constants::nFv0Channels; iPmt++) {
  for (int iPmt = 0; iPmt < 16; iPmt++) {
    // run the CFD: this updates the start index for the next BC in case the CFD dead time ends in the next BC
    double cfdWithOffset = SimulateTimeCfd(mCfdStartIndex[iPmt], mLastBCCache.mPmtChargeVsTime[iPmt], bc.mPmtChargeVsTime[iPmt]);
    double cfdZero = cfdWithOffset - FD2DigParam::Instance().avgCfdTimeForMip;

    // Conditions to sum charge are: all participating channels must have time within +/- 2.5 ns, AND
    //   at least one channel must follow more strict conditions (see below)
    if (cfdZero < -FD2DigParam::Instance().cfdCheckWindow || cfdZero > FD2DigParam::Instance().cfdCheckWindow) {
      continue;
    }

    int iTotalCharge = std::lround(IntegrateCharge(bc.mPmtChargeVsTime[iPmt]) * DP::INV_CHARGE_PER_ADC); // convert Coulomb to adc;

    uint8_t channelBits = FD2DigParam::Instance().defaultChainQtc;
    if (std::rand() % 2) {
       o2::fv0::ChannelData::setFlag(o2::fv0::ChannelData::kNumberADC, channelBits);
    }
    if (iTotalCharge > (FD2DigParam::Instance().maxCountInAdc) && FD2DigParam::Instance().useMaxChInAdc) {
      iTotalCharge = FD2DigParam::Instance().maxCountInAdc; // max adc channel for one PMT
      o2::fv0::ChannelData::setFlag(o2::fv0::ChannelData::kIsAmpHigh, channelBits);
    }

    if (iTotalCharge < FD2DigParam::Instance().getCFDTrshInAdc()) {
      continue;
    }

    int iCfdZero = std::lround(cfdZero * DP::INV_TIME_PER_TDCCHANNEL);
    digitsCh.emplace_back(iPmt, iCfdZero, iTotalCharge, channelBits);
    ++nTotFiredCells;

    int triggerGate = FD2DigParam::Instance().mTime_trg_gate;
    if (std::abs(iCfdZero) < triggerGate) {
      ++nTrgFiredCells;
      //---trigger---
      totalChargeAllRing += iTotalCharge;
      avgTime += iCfdZero;
      if (iPmt < 24) {
        nSignalInner++;
        totalChargeInnerRing += iTotalCharge;
      } else {
        nSignalOuter++;
        totalChargeOuterRing += iTotalCharge;
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
    avgTime = o2::fit::Triggers::DEFAULT_TIME;
  }
  /// Triggers for FD2
  bool isA, isNchannels, isAIn, isAOut, isTotalCharge;
  isA = nTrgFiredCells > 0;
  isNchannels = nTrgFiredCells > FD2DigParam::Instance().NchannelsLevel;
  // isAIn = nSignalInner > FD2DigParam::Instance().NchannelsLevel;  // ring 1,2 and 3
  isAIn = 0.125 * totalChargeInnerRing > 2 * FD2DigParam::Instance().InnerChargeLevel; // ring 1,2 and 3
  // isAOut = nSignalOuter > FD2DigParam::Instance().NchannelsLevel; // ring 4 and 5
  isAOut = 0.125 * totalChargeOuterRing > 2 * FD2DigParam::Instance().OuterChargeLevel; // ring 4 and 5
  isTotalCharge = 0.125 * totalChargeAllRing > 2 * FD2DigParam::Instance().ChargeLevel;

  o2::fit::Triggers triggers;
  const int unusedCharge = o2::fit::Triggers::DEFAULT_AMP;
  const int unusedTime = o2::fit::Triggers::DEFAULT_TIME;
  const int unusedZero = o2::fit::Triggers::DEFAULT_ZERO;
  const bool unusedBitsInSim = false; // bits related to laser and data validity
  const bool bitDataIsValid = true;
  triggers.setTriggers(isA, isAIn, isAOut, isTotalCharge, isNchannels, nTrgFiredCells, (int8_t)unusedZero,
                       (int32_t)(0.125 * totalChargeAllRing), (int32_t)unusedCharge, (int16_t)avgTime, (int16_t)unusedTime, unusedBitsInSim, unusedBitsInSim, bitDataIsValid);
  digitsBC.emplace_back(first, nTotFiredCells, bc, triggers, mEventId - 1);
  digitsTrig.emplace_back(bc, isA, isAIn, isAOut, isTotalCharge, isNchannels);
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
  const Float_t epsilon = 0.0001f;
  const Float_t p = FD2DigParam::Instance().lightYield * FD2DigParam::Instance().photoCathodeEfficiency;
  if ((fabs(1.0f - p) < epsilon) || nPhot == 0) {
    return nPhot;
  }
  const Int_t n = Int_t(nPhot < 100
                          ? gRandom->Binomial(nPhot, p)
                          : gRandom->Gaus((p * nPhot) + 0.5, TMath::Sqrt(p * (1. - p) * nPhot)));
  return n;
}
//---------------------------------------------------------------------------
Float_t Digitizer::IntegrateCharge(const ChannelDigitF& pulse) const
{
  int const chargeIntMin = FD2DigParam::Instance().isIntegrateFull ? 0 : (FD2DigParam::Instance().avgCfdTimeForMip - 6.0) / mBinSize;                // Charge integration offset (cfd mean time - 6 ns)
  int const chargeIntMax = FD2DigParam::Instance().isIntegrateFull ? mNTimeBinsPerBC : (FD2DigParam::Instance().avgCfdTimeForMip + 14.0) / mBinSize; // Charge integration offset (cfd mean time + 14 ns)
  if (chargeIntMin < 0 || chargeIntMin > mNTimeBinsPerBC || chargeIntMax > mNTimeBinsPerBC) {
    LOG(fatal) << "invalid indicess: chargeInMin=" << chargeIntMin << " chargeIntMax=" << chargeIntMax;
  }
  Float_t totalCharge = 0.0f;
  for (int iTimeBin = chargeIntMin; iTimeBin < chargeIntMax; iTimeBin++) {
    totalCharge += pulse[iTimeBin];
  }
  return totalCharge;
}
//---------------------------------------------------------------------------
Float_t Digitizer::SimulateTimeCfd(int& startIndex, const ChannelDigitF& pulseLast, const ChannelDigitF& pulse) const
{
  Float_t timeCfd = -1024.0f;

  if (pulse.empty()) {
    startIndex = 0;
    return timeCfd;
  }

  Float_t const cfdThrInCoulomb = FD2DigParam::Instance().mCFD_trsh * 1e-3 / 50 * mBinSize * 1e-9; // convert mV into Coulomb assuming 50 Ohm

  Int_t const binShift = TMath::Nint(FD2DigParam::Instance().timeShiftCfd / mBinSize);
  Float_t sigPrev = 5 * pulseLast[mNTimeBinsPerBC - binShift - 1] - pulseLast[mNTimeBinsPerBC - 1]; //  CFD output from the last bin of the last BC
  for (Int_t iTimeBin = 0; iTimeBin < mNTimeBinsPerBC; ++iTimeBin) {
    Float_t const sigCurrent = 5.0f * (iTimeBin >= binShift ? pulse[iTimeBin - binShift] : pulseLast[mNTimeBinsPerBC - binShift + iTimeBin]) - pulse[iTimeBin];
    if (iTimeBin >= startIndex && std::abs(pulse[iTimeBin]) > cfdThrInCoulomb) { // enable
      if (sigPrev < 0.0f && sigCurrent >= 0.0f) {                                // test for zero-crossing
        timeCfd = Float_t(iTimeBin) * mBinSize;
        startIndex = iTimeBin + std::lround(FD2DigParam::Instance().mCfdDeadTime / mBinSize); // update startIndex (CFD dead time)
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

float Digitizer::getDistFromCellCenter(UInt_t cellId, double hitx, double hity)
{
 // Geometry* geo = Geometry::instance();


  // Parametrize the line (ax+by+c=0) that crosses the detector center and the cell's middle point
  // Point3Dsimple* pCell = &geo->getCellCenter(cellId);

  Point3Dsimple* pCell = mCellCenters.at(cellId);	

  // float x0, y0, z0;
  // geo->getGlobalPosition(x0, y0, z0);
  float x0 = 0, y0 = 9, z0 = cellId < 8 ? Constants::zScintA : Constants::zScintC;
  double a = -(y0 - pCell->y) / (x0 - pCell->x);
  double b = 1;
  double c = -(y0 - a * x0);
  // Return the distance from hit to this line
  return (a * hitx + b * hity + c) / TMath::Sqrt(a * a + b * b);
}

float Digitizer::getSignalFraction(float distanceFromXc, bool isFirstChannel)
{
  float const fraction = sigmoidPmtRing5(distanceFromXc);
  if (distanceFromXc > 0) {
    return isFirstChannel ? fraction : (1. - fraction);
  } else {
    return isFirstChannel ? (1. - fraction) : fraction;
  }
}

//_____________________________________________________________________________
o2::fd2::Digitizer::BCCache& Digitizer::setBCCache(const o2::InteractionRecord& ir)
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
o2::fd2::Digitizer::BCCache* Digitizer::getBCCache(const o2::InteractionRecord& ir)
{
  // get pointer on existing cache
  for (auto cb = mCache.begin(); cb != mCache.end(); cb++) {
    if ((*cb) == ir) {
      return &(*cb);
    }
  }
  return nullptr;
}

bool Digitizer::isRing5(int detID)
{
  if (detID > 31) {
    return true;
  } else {
    return false;
  }
}

//_____________________________________________________________________________
void Digitizer::InitCellCenters()
{
  auto& baseParam = o2::fd2::FD2BaseParam::Instance();

  for (int cellId = 0; cellId < 16; cellId ++) {
    float rmin = cellId < 8 ?  baseParam.rminA : baseParam.rminC;
    float rmax = cellId < 8 ?  baseParam.rmaxA : baseParam.rmaxC;
    float rc = 0.5 * (rmin + rmax);
     
    int ic = cellId % 8;
    float theta = (2 * ic +1) * TMath::Pi() / 8;
    float x = rc * TMath::Cos(theta);
    float y = rc * TMath::Sin(theta);

    float z = cellId < 8 ? Constants::zScintA : Constants::zScintC;

    Point3Dsimple* p;
    p->SetCoordinates(x, y, z);
    
    mCellCenters[cellId] = p;
  }
}

O2ParamImpl(FD2DigParam);
