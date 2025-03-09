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

/// \file Detector.cxx
/// \brief Implementation of the Detector class

#include "ITSMFTSimulation/Hit.h"
#include "FDSimulation/Detector.h"
#include "FDBase/GeometryTGeo.h"
#include "FDBase/FDBaseParam.h"
#include "FDBase/Constants.h"

#include "DetectorsBase/Stack.h"
#include "SimulationDataFormat/TrackReference.h"
#include "Field/MagneticField.h"

// FairRoot includes
#include "FairDetector.h"
#include <fairlogger/Logger.h>
#include "FairRootManager.h"
#include "FairRun.h"
#include "FairRuntimeDb.h"
#include "FairVolume.h"
#include "FairRootManager.h"

#include "TVirtualMC.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <TGeoTube.h>
#include <TGeoVolume.h>
#include <TGeoCompositeShape.h>
#include <TGeoMedium.h>
#include <TGeoCone.h>
#include <TGeoManager.h>
#include "TRandom.h"

class FairModule;

class TGeoMedium;

using namespace o2::fd;
using o2::itsmft::Hit;

Detector::Detector(bool active)
  : o2::base::DetImpl<Detector>("FD", true),
    mHits(o2::utils::createSimVector<o2::itsmft::Hit>()),
    mGeometryTGeo(nullptr),
    mTrackData()
{
  mNumberOfRingsC = Constants::nringsC;
  mNumberOfSectors = Constants::nsect;

  mEtaMinA = Constants::etaMin;
  mEtaMaxA = Constants::etaMax;
  mEtaMinC = -Constants::etaMax;
  mEtaMaxC = -Constants::etaMin;

  auto& baseParam = FDBaseParam::Instance();

  if (baseParam.withMG) {
    mNumberOfRingsA = Constants::nringsA_withMG;
    mEtaMinA = Constants::etaMinA_withMG;
  } else {
    mNumberOfRingsA = Constants::nringsA;
    mEtaMinA = Constants::etaMin;
  }

  mDzScint = baseParam.dzscint / 2;
  mZmodA = baseParam.zmodA;
  mZmodC = baseParam.zmodC;

  for (int i = 0; i <= mNumberOfRingsA + 1; i++) {
    float eta = mEtaMaxA - i * (mEtaMaxA - mEtaMinA) / mNumberOfRingsA;
    float r = ringRadius(mZmodA, eta);
    mRingRadiiA.emplace_back(r);
  }

  for (int i = 0; i <= mNumberOfRingsC + 1; i++) {
    float eta = mEtaMaxC - i * (mEtaMaxC - mEtaMinC) / mNumberOfRingsC;
    float r = ringRadius(mZmodC, eta);
    mRingRadiiC.emplace_back(r);
  }
}

Detector::Detector(const Detector& rhs)
  : o2::base::DetImpl<Detector>(rhs),
    mTrackData(),
    mHits(o2::utils::createSimVector<o2::itsmft::Hit>())
{
}

Detector& Detector::operator=(const Detector& rhs)
{

  if (this == &rhs) {
    return *this;
  }
  // base class assignment
  base::Detector::operator=(rhs);
  mTrackData = rhs.mTrackData;

  mHits = nullptr;
  return *this;
}

Detector::~Detector()
{

  if (mHits) {
    o2::utils::freeSimVector(mHits);
  }
}

void Detector::InitializeO2Detector()
{
  LOG(info) << "Initialize FD detector";
  mGeometryTGeo = GeometryTGeo::Instance();
  defineSensitiveVolumes();
}

bool Detector::ProcessHits(FairVolume* vol)
{
  // This method is called from the MC stepping
  if (!(fMC->TrackCharge())) {
    return kFALSE;
  }

  auto stack = (o2::data::Stack*)fMC->GetStack();

  // int cellId = vol->getVolumeId();

  // Check track status to define when hit is started and when it is stopped
  bool startHit = false, stopHit = false;
  unsigned char status = 0;
  if (fMC->IsTrackEntering()) {
    status |= Hit::kTrackEntering;
  }
  if (fMC->IsTrackInside()) {
    status |= Hit::kTrackInside;
  }
  if (fMC->IsTrackExiting()) {
    status |= Hit::kTrackExiting;
  }
  if (fMC->IsTrackOut()) {
    status |= Hit::kTrackOut;
  }
  if (fMC->IsTrackStop()) {
    status |= Hit::kTrackStopped;
  }
  if (fMC->IsTrackAlive()) {
    status |= Hit::kTrackAlive;
  }

  // track is entering or created in the volume
  if ((status & Hit::kTrackEntering) || (status & Hit::kTrackInside && !mTrackData.mHitStarted)) {
    startHit = true;
  } else if ((status & (Hit::kTrackExiting | Hit::kTrackOut | Hit::kTrackStopped))) {
    stopHit = true;
  }

  // increment energy loss at all steps except entrance
  if (!startHit) {
    mTrackData.mEnergyLoss += fMC->Edep();
  }
  if (!(startHit | stopHit)) {
    return kFALSE; // do noting
  }

  if (startHit) {
    mTrackData.mHitStarted = true;
    mTrackData.mEnergyLoss = 0.;
    fMC->TrackMomentum(mTrackData.mMomentumStart);
    fMC->TrackPosition(mTrackData.mPositionStart);
    mTrackData.mTrkStatusStart = true;
  }

  if (stopHit) {
    TLorentzVector positionStop;
    fMC->TrackPosition(positionStop);
    int trackId = stack->GetCurrentTrackNumber();
    unsigned int chId = getChannelId(mTrackData.mPositionStart.Vect());

    Hit* p = addHit(trackId, chId /*cellId*/, mTrackData.mPositionStart.Vect(), positionStop.Vect(),
                    mTrackData.mMomentumStart.Vect(), mTrackData.mMomentumStart.E(),
                    positionStop.T(), mTrackData.mEnergyLoss, mTrackData.mTrkStatusStart,
                    status);
    stack->addHit(GetDetId());
  } else {
    return false; // do nothing more
  }
  return true;
}

o2::itsmft::Hit* Detector::addHit(int trackId, int cellId,
                                  const TVector3& startPos,
                                  const TVector3& endPos,
                                  const TVector3& startMom,
                                  double startE,
                                  double endTime,
                                  double eLoss,
                                  unsigned int startStatus,
                                  unsigned int endStatus)
{
  mHits->emplace_back(trackId, cellId, startPos,
                      endPos, startMom, startE, endTime, eLoss, startStatus, endStatus);
  return &(mHits->back());
}

void Detector::ConstructGeometry()
{
  createMaterials();
  buildModules();
}

void Detector::EndOfEvent()
{
  Reset();
}

void Detector::Register()
{
  // This will create a branch in the output tree called Hit, setting the last
  // parameter to kFALSE means that this collection will not be written to the file,
  // it will exist only during the simulation

  if (FairRootManager::Instance()) {
    FairRootManager::Instance()->RegisterAny(addNameTo("Hit").data(), mHits, kTRUE);
  }
}

void Detector::Reset()
{
  if (!o2::utils::ShmManager::Instance().isOperational()) {
    mHits->clear();
  }
}

void Detector::createMaterials()
{

  float density, as[11], zs[11], ws[11];
  double radLength, absLength, a_ad, z_ad;
  int id;

  // EJ-204 scintillator, based on polyvinyltoluene
  const int nScint = 2;
  float aScint[nScint] = {1.00784, 12.0107};
  float zScint[nScint] = {1, 6};
  float wScint[nScint] = {0.07085, 0.92915}; // based on EJ-204 datasheet: n_atoms/cm3
  const float dScint = 1.023;

  int matId = 0;                  // tmp material id number
  const int unsens = 0, sens = 1; // sensitive or unsensitive medium
                                  //
  int fieldType = 3;              // Field type
  float maxField = 5.0;           // Field max.

  float tmaxfd = -10.0; // max deflection angle due to magnetic field in one step
  float stemax = 0.1;   // max step allowed [cm]
  float deemax = 1.0;   // maximum fractional energy loss in one step 0<deemax<=1
  float epsil = 0.03;   // tracking precision [cm]
  float stmin = -0.001; // minimum step due to continuous processes [cm] (negative value: choose it automatically)

  LOG(info) << "FD: CreateMaterials(): fieldType " << fieldType << ", maxField " << maxField;

  o2::base::Detector::Mixture(++matId, "Scintillator", aScint, zScint, dScint, nScint, wScint);
  o2::base::Detector::Medium(Scintillator, "Scintillator", matId, unsens, fieldType, maxField,
                             tmaxfd, stemax, deemax, epsil, stmin);
}

void Detector::buildModules()
{
  LOGP(info, "Creating FD geometry");

  TGeoVolume* vCave = gGeoManager->GetVolume("cave");

  if (!vCave) {
    LOG(fatal) << "Could not find the top volume (cave)!";
  }

  TGeoVolumeAssembly* vFDA = buildModuleA();
  TGeoVolumeAssembly* vFDC = buildModuleC();

  vCave->AddNode(vFDA, 1, new TGeoTranslation(0., 0., mZmodA /* - mDzScint/2.*/));
  vCave->AddNode(vFDC, 1, new TGeoTranslation(0., 0., mZmodC /* + mDzScint/2.*/));
}

TGeoVolumeAssembly* Detector::buildModuleA()
{
  TGeoVolumeAssembly* mod = new TGeoVolumeAssembly("FDA");

  const TGeoMedium* medium = gGeoManager->GetMedium("FD_Scintillator");

  float dphiDeg = 360. / mNumberOfSectors;

  for (int ir = 0; ir < mNumberOfRingsA; ir++) {
    std::string rName = "fd_ring" + std::to_string(ir + 1);
    TGeoVolumeAssembly* ring = new TGeoVolumeAssembly(rName.c_str());
    float rmin = mRingRadiiA[ir];
    float rmax = mRingRadiiA[ir + 1];
    LOG(info) << "ring" << ir << ": from " << rmin << " to " << rmax;
    for (int ic = 0; ic < mNumberOfSectors; ic++) {
      int cellId = ic + mNumberOfSectors * ir;
      std::string nodeName = "fd_node" + std::to_string(cellId);
      float phimin = dphiDeg * ic;
      float phimax = dphiDeg * (ic + 1);
      auto tbs = new TGeoTubeSeg("tbs", rmin, rmax, mDzScint, phimin, phimax);
      auto nod = new TGeoVolume(nodeName.c_str(), tbs, medium);
      nod->SetLineColor(kRed);
      ring->AddNode(nod, cellId);
    }
    mod->AddNode(ring, 1);
  }

  return mod;
}

TGeoVolumeAssembly* Detector::buildModuleC()
{
  TGeoVolumeAssembly* mod = new TGeoVolumeAssembly("FDC");

  const TGeoMedium* medium = gGeoManager->GetMedium("FD_Scintillator");

  float dphiDeg = 360. / mNumberOfSectors;

  for (int ir = 0; ir < mNumberOfRingsC; ir++) {
    std::string rName = "fd_ring" + std::to_string(ir + 1 + mNumberOfRingsA);
    TGeoVolumeAssembly* ring = new TGeoVolumeAssembly(rName.c_str());
    float rmin = mRingRadiiC[ir];
    float rmax = mRingRadiiC[ir + 1];
    LOG(info) << "ring" << ir + mNumberOfRingsA << ": from " << rmin << " to " << rmax;
    for (int ic = 0; ic < mNumberOfSectors; ic++) {
      int cellId = ic + mNumberOfSectors * (ir + mNumberOfRingsA);
      std::string nodeName = "fd_node" + std::to_string(cellId);
      float phimin = dphiDeg * ic;
      float phimax = dphiDeg * (ic + 1);
      auto tbs = new TGeoTubeSeg("tbs", rmin, rmax, mDzScint, phimin, phimax);
      auto nod = new TGeoVolume(nodeName.c_str(), tbs, medium);
      nod->SetLineColor(kBlue);
      ring->AddNode(nod, cellId);
    }
    mod->AddNode(ring, 1);
  }

  return mod;
}

void Detector::defineSensitiveVolumes()
{
  LOG(info) << "Adding FD Sentitive Volumes";

  TGeoVolume* v;
  TString volumeName;

  int nCellA = mNumberOfRingsA * mNumberOfSectors;
  int nCellC = mNumberOfRingsC * mNumberOfSectors;

  LOG(info) << "number of A rings = " << mNumberOfRingsA << " number of cells = " << nCellA;

  for (int iv = 0; iv < nCellA + nCellC; iv++) {
    volumeName = "fd_node" + std::to_string(iv);
    v = gGeoManager->GetVolume(volumeName);
    LOG(info) << "Adding FD Sensitive Volume => " << v->GetName();
    AddSensitiveVolume(v);
  }
}

unsigned int Detector::getChannelId(TVector3 vec)
{
  float phi = vec.Phi();
  if (phi < 0) {
    phi += TMath::TwoPi();
  }

  float r = vec.Perp();
  float z = vec.Z();

  int isect = int(phi / (TMath::TwoPi() / mNumberOfSectors));

  std::vector<float> rd = z > 0 ? mRingRadiiA : mRingRadiiC;
  int noff = z > 0 ? 0 : mNumberOfRingsA * mNumberOfSectors;

  int ir = 0;

  for (int i = 1; i < rd.size(); i++) {
    if (r < rd[i]) {
      break;
    } else {
      ir++;
    }
  }

  return ir * mNumberOfSectors + isect + noff;
}

float Detector::ringRadius(float z, float eta)
{
  return z * TMath::Tan(2 * TMath::ATan(TMath::Exp(-eta)));
}

ClassImp(o2::fd::Detector);
