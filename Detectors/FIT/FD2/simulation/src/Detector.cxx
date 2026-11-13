// Copyright 2019-2026 CERN and copyright holders of ALICE O2.
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

#include "FD2Simulation/Detector.h"
#include "FD2Base/GeometryTGeo.h"
#include "FD2Base/FD2BaseParam.h"
#include "FD2Base/Constants.h"

#include "DetectorsBase/Stack.h"
#include "SimulationDataFormat/TrackReference.h"
#include "Field/MagneticField.h"

// FairRoot includes
#include "FairDetector.h"
#include <fairlogger/Logger.h>
#include "FairRootManager.h"
#include "FairRun.h"
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
#include <cmath>

class FairModule;

class TGeoMedium;

using namespace o2::fd2;
using o2::fd3::Hit;

Detector::Detector(bool active)
  : o2::base::DetImpl<Detector>("FD2", true),
    mHits(o2::utils::createSimVector<o2::fd3::Hit>()),
    mGeometryTGeo(nullptr),
    mTrackData()
{
}

Detector::Detector(const Detector& rhs)
  : o2::base::DetImpl<Detector>(rhs),
    mTrackData(),
    mHits(o2::utils::createSimVector<o2::fd3::Hit>()),
    mChannelCounter(0),
    mChannelId{}
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
  LOG(info) << "Initialize Forward Detector";
  mGeometryTGeo = GeometryTGeo::Instance();
  defineSensitiveVolumes();
  definePassiveVolumes();
}

bool Detector::ProcessHits(FairVolume* vol)
{
  // This method is called from the MC stepping
  if (!(fMC->TrackCharge())) {
    return kFALSE;
  }

  int copy;
  int volId = fMC->CurrentVolID(copy);
  int detId = mChannelId[volId];

  auto stack = (o2::data::Stack*)fMC->GetStack();

  // Check track status to define when hit is started and when it is stopped
  int particlePdg = fMC->TrackPid();
  bool startHit = false, stopHit = false;
  if ((fMC->IsTrackEntering()) || (fMC->IsTrackInside() && !mTrackData.mHitStarted)) {
    startHit = true;
  } else if ((fMC->IsTrackExiting() || fMC->IsTrackOut() || fMC->IsTrackStop())) {
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

    math_utils::Point3D<float> posStart(mTrackData.mPositionStart.X(), mTrackData.mPositionStart.Y(), mTrackData.mPositionStart.Z());
    math_utils::Point3D<float> posStop(positionStop.X(), positionStop.Y(), positionStop.Z());
    math_utils::Vector3D<float> momStart(mTrackData.mMomentumStart.Px(), mTrackData.mMomentumStart.Py(), mTrackData.mMomentumStart.Pz());

    Hit* p = addHit(trackId, detId, posStart, posStop,
                    momStart, mTrackData.mMomentumStart.E(),
                    positionStop.T(), mTrackData.mEnergyLoss, particlePdg);
    stack->addHit(GetDetId());
    mTrackData.mHitStarted = false;
  } else {
    return false; // do nothing more
  }
  return true;
}

o2::fd3::Hit* Detector::addHit(int trackId, unsigned int detId,
                               const math_utils::Point3D<float>& startPos,
                               const math_utils::Point3D<float>& endPos,
                               const math_utils::Vector3D<float>& startMom,
                               double startE,
                               double endTime,
                               double eLoss,
                               int particlePdg)
{
  mHits->emplace_back(trackId, detId, startPos,
                      endPos, startMom, startE, endTime, eLoss, particlePdg);
  return &(mHits->back());
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

void Detector::ConstructGeometry()
{
  createMaterials();
  buildModules();
}

void Detector::createMaterials()
{
  int matId = 0;                  // tmp material id number
  const int unsens = 0, sens = 1; // sensitive or unsensitive medium
  int fieldType;
  float maxField;

  // EJ-204 scintillator, based on polyvinyltoluene
  float aScint[2] = {1.00784, 12.0107};
  float zScint[2] = {1, 6};
  float wScint[2] = {0.07085, 0.92915}; // based on EJ-204 datasheet: n_atoms/cm3
  const float dScint = 1.023;

  // Radiator  glass SiO2
  float aglass[2] = {28.0855, 15.9994};
  float zglass[2] = {14., 8.};
  float wglass[2] = {1., 2.};
  const float dglass = 2.2;

  // PMMA plastic mixture: (C5O2H8)n, same for plastic fiber support and for the fiber core
  //   Fiber cladding is different, but it comprises only 3% of the fiber volume, so it is not included
  const Int_t nPlast = 3;
  Float_t aPlast[nPlast] = {1.00784, 12.0107, 15.999};
  Float_t zPlast[nPlast] = {1, 6, 8};
  Float_t wPlast[nPlast] = {0.08054, 0.59985, 0.31961};
  const Float_t dPlast = 1.18;

  o2::base::Detector::initFieldTrackingParams(fieldType, maxField);
  // TODO: Comment out two lines below once tested that the above function assigns field type and max correctly
  fieldType = 3;  // Field type
  maxField = 5.0; // Field max.

  float tmaxfd = -10.0;  // max deflection angle due to magnetic field in one step
  float stepmax = 0.1;    // max step allowed [cm]
  float deemax = 1.0;     // maximum fractional energy loss in one step 0<deemax<=1
  float epsil = 0.03;     // tracking precision [cm]
  float stepmin = -0.001; // minimum step due to continuous processes [cm] (negative value: choose it automatically)

  o2::base::Detector::Mixture(++matId, "Scintillator", aScint, zScint, dScint, 2, wScint);
  o2::base::Detector::Medium(Scintillator, "Scintillator", matId, sens, fieldType, maxField,
                             tmaxfd, stepmax, deemax, epsil, stepmin);

  // Cherenkov radiator glass, modify fieldType and maxField parameters
  fieldType = 2;
  maxField = 10.;

  o2::base::Detector::Mixture(++matId, "Radiator optical glass", aglass, zglass, dglass, -2, wglass);
  o2::base::Detector::Medium(RadiatorOpticalGlass, "OpticalGlass", matId, sens, fieldType, maxField,
                             10., .01, .1, .003, .01);

  o2::base::Detector::Mixture(++matId, "Fiber", aPlast, zPlast, /*dFiber*/ dPlast, nPlast, wPlast);
  o2::base::Detector::Medium(Fiber, "Fiber", matId, unsens, fieldType, maxField,
                               tmaxfd, stepmax, deemax, epsil, stepmin);
}

void Detector::buildModules()
{
  LOGP(info, "Creating FD2 geometry");

  auto& baseParam = FD2BaseParam::Instance();

  auto vCaveRB24 = (TGeoVolume*)gGeoManager->GetVolume("caveRB24");
  auto vCave     = (TGeoVolume*)gGeoManager->GetVolume("cave");

  mChannelCounter = 0;

  TGeoVolumeAssembly *vFD2_ScintA, *vFD2_ScintC;

   float rminA = baseParam.rminA, rmaxA = baseParam.rmaxA,
	 rminC = baseParam.rminC, rmaxC = baseParam.rmaxC;

  //vFD2_ScintA = buildModuleScint(Constants::zScintA, Constants::etaMinA, Constants::etaMaxA);
  vFD2_ScintA = buildModuleScint(rminA, rmaxA);
  vFD2_ScintC = buildModuleScint(rminC, rmaxC);

  vFD2_ScintA->SetName("FD2_ScintA");
  vFD2_ScintC->SetName("FD2_ScintC");

 // vFD2_ScintC = buildModuleScint(Constants::zScintC, Constants::etaMaxC, Constants::etaMinC);
  vCaveRB24->AddNode(vFD2_ScintA, 1, new TGeoTranslation(0.f, 0.f, Constants::zScintA  - Constants::offsetRB24)); 
  vCave->AddNode(vFD2_ScintC, 1, new TGeoTranslation(0.f, 0.f, Constants::zScintC));

 // float rA = 0.5 *(getRingSize(Constants::zScintA, Constants::etaMaxA) + getRingSize(Constants::zScintA, Constants::etaMinA));
  float rA = 0.5 * (rmaxA + rminA);
  for (int i = 0; i < Constants::nSectScint; i ++) {
    auto fb = buildFiberBundle();
    std::string fbname = "bundle" + std::to_string(i) + "_A";
    fb->SetName(fbname.c_str());
    float theta = (2 * i + 1) * TMath::Pi() / Constants::nSectScint ;
    vCaveRB24->AddNode(fb, 1, new TGeoTranslation(rA * TMath::Cos(theta), rA * TMath::Sin(theta), Constants::zScintA - Constants::offsetRB24 + Constants::dzScint / 2 + Constants::lenFiber / 2));
  }

  // float rC = 0.5 *(getRingSize(Constants::zScintC, Constants::etaMaxC) + getRingSize(Constants::zScintC, Constants::etaMinC));
  float rC = 0.5 * (rmaxC + rminC);
  for (int i = 0; i < Constants::nSectScint; i ++) {
    auto fb = buildFiberBundle();
    std::string fbname = "bundle" + std::to_string(i) + "_C";
    fb->SetName(fbname.c_str());
    float theta = (2 * i + 1) * TMath::Pi() / Constants::nSectScint ;
    vCave->AddNode(fb, 1, new TGeoTranslation(rC * TMath::Cos(theta), rC * TMath::Sin(theta), Constants::zScintC - Constants::dzScint / 2 - Constants::lenFiber /2));
  }
}

TGeoVolumeAssembly* Detector::buildModuleScint(float zmod, float etaMin, float etaMax)
{
  auto mod = new TGeoVolumeAssembly("");

  const TGeoMedium* medium = gGeoManager->GetMedium("FD2_Scintillator");

  int nrings = Constants::nRingsScint, nsect = Constants::nSectScint;
  float dz = Constants::dzScint;

  float r1 = getRingSize(zmod, etaMax), r2 = getRingSize(zmod, etaMin);
  if (r1 >= r2) {
    LOG(fatal) << "rmin(" << r1 << ") should be less than rmax(" << r2 << "). Please check your input!";
  }
  LOG(info) << "Scintillator ring: from " << r1 << " to " << r2;
  for (int ic = 0; ic < nsect; ic++) {
    int cellId = mChannelCounter++;
    std::string nodeName = "fd2_node" + std::to_string(cellId);
    float phi1 = ic * 360. / nsect;
    float phi2 = (ic + 1) * 360. / nsect;
    auto tbs = new TGeoTubeSeg(r1, r2, dz / 2, phi1, phi2);
    auto node = new TGeoVolume(nodeName.c_str(), tbs, medium);
    if (cellId < nrings * nsect) {
      node->SetLineColor(ic % 2 == 0 ? kRed : kRed - 7);
    } else {
      node->SetLineColor(ic % 2 == 0 ? kBlue : kBlue - 7);
    }
    mod->AddNode(node, 1);
  }

  return mod;
}

TGeoVolumeAssembly* Detector::buildModuleScint(float rmin, float rmax)
{
  if (rmin >= rmax) {
    LOG(fatal) << "rmin(" << rmin << ") should be less than rmax(" << rmax << "). Please check your input!";
  }

  auto mod = new TGeoVolumeAssembly("");

  const TGeoMedium* medium = gGeoManager->GetMedium("FD2_Scintillator");

  int   nsect = Constants::nSectScint;
  float dz = Constants::dzScint;

  LOG(info) << "Scintillator ring: from " << rmin << " to " << rmax;
  for (int ic = 0; ic < nsect; ic++) {
    int cellId = mChannelCounter++;
    std::string nodeName = "fd2_node" + std::to_string(cellId);
    float phi1 = ic * 360. / nsect;
    float phi2 = (ic + 1) * 360. / nsect;
    auto tbs = new TGeoTubeSeg(rmin, rmax, dz / 2, phi1, phi2);
    auto node = new TGeoVolume(nodeName.c_str(), tbs, medium);
    if (cellId <  nsect) {
      node->SetLineColor(ic % 2 == 0 ? kRed : kRed - 7);
    } else {
      node->SetLineColor(ic % 2 == 0 ? kBlue : kBlue - 7);
    }
    mod->AddNode(node, 1);
  }

  return mod;
}

TGeoVolumeAssembly* Detector::buildFiberBundle()
{
  auto bundle = new TGeoVolumeAssembly("");
  
  const TGeoMedium* medium = gGeoManager->GetMedium("FD2_Fiber");

  std::vector<std::pair<float, float>> xyFiber;
  float dxy = 0.2; // fiber spacing

  for (int i = 0; i < 196; i ++) {
     int ir = i / 14,  ic = i % 14;
     float x = - dxy * 7 + dxy * ir;
     float y =   dxy * 7 - dxy * ic;
     xyFiber.emplace_back(x, y);
  }
//  xyFiber.emplace_back(7.5 * dxy, 7.5 * dxy);
//  xyFiber.emplace_back(7.5 * dxy, -7.5 * dxy);
//  xyFiber.emplace_back(-7.5 * dxy, -7.5 * dxy);
//  xyFiber.emplace_back(-7.5 * dxy, 7.5 * dxy);

  int ii = 0;
  for (auto& xyf : xyFiber) {
    auto tube = new TGeoTube(Constants::rMinFiber, Constants::rMaxFiber, Constants::lenFiber / 2);
    std::string nodeName = "fiber" + std::to_string(ii);
    auto node = new TGeoVolume(nodeName.c_str(), tube, medium);
    node->SetLineColor(kOrange + 3);
    bundle->AddNode(node, 1, new TGeoTranslation(xyf.first, xyf.second, 0.f));
    ii ++;
  }

  return bundle;
}

void Detector::defineSensitiveVolumes()
{
  LOG(info) << "Adding FD2 sentitive volumes...";

  mChannelId = {};

  auto FD2_ScintA =  (TGeoVolumeAssembly*)gGeoManager->GetVolume("FD2_ScintA");
  auto FD2_ScintC =  (TGeoVolumeAssembly*)gGeoManager->GetVolume("FD2_ScintC");

  int ivol = 0;

  for (int i = 0; i < FD2_ScintA->GetNdaughters(); i++) {
     auto node = (TGeoNode*)FD2_ScintA->GetNode(i);
     auto v = (TGeoVolume*)node->GetVolume();
     int volId = registerSensitiveVolumeAndGetVolID(v);
     LOG(info) << "volId = " << volId << " ivol = " << ivol; 
     mChannelId[volId] = ivol++;
  }

  for (int i = 0; i < FD2_ScintC->GetNdaughters(); i++) {
     auto node = (TGeoNode*)FD2_ScintC->GetNode(i);
     auto v = (TGeoVolume*)node->GetVolume();
     int volId = registerSensitiveVolumeAndGetVolID(v);
     LOG(info) << "volId = " << volId << " ivol = " << ivol; 
     mChannelId[volId] = ivol++;
  }

//  for (int ivol = 0; ivol < mChannelCounter; ivol++) {
//    std::string volumeName = "fd2_node" + std::to_string(ivol);
//    auto v = (TGeoVolume*)gGeoManager->GetVolume(volumeName.c_str());
//    if (!v) {
//      continue;
//    }
//    int volId = registerSensitiveVolumeAndGetVolID(v);
//    mChannelId[volId] = ivol;
//  }
  LOG(info) << "Done";
}

void Detector::definePassiveVolumes()
{
  // To be added later
}

float Detector::getRingSize(float z, float eta)
{
  return z * TMath::Tan(2 * TMath::ATan(TMath::Exp(-eta)));
}

ClassImp(o2::fd2::Detector);
