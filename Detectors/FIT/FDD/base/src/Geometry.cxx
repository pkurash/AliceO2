// }
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
#include <iomanip>
#include "FDDBase/Geometry.h"
#include "FDDBase/Constants.h"
#include "FDDBase/FDDBaseParam.h"

#include <TGeoManager.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoCompositeShape.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoMatrix.h>
#include "Framework/Logger.h"

#include "TMath.h"

ClassImp(o2::fdd::Geometry);

using namespace o2::fdd;
//_____________________________________________________________________________
Geometry::Geometry(EGeoType initType)
{
  mGeometryType = initType;
  buildGeometry();
}
//_____________________________________________________________________________
Geometry::Geometry(const Geometry& geom)
{
  this->mGeometryType = geom.mGeometryType;
}
//_____________________________________________________________________________
void Geometry::buildGeometry()
{

  LOG(info) << "Geometry::buildGeometry()::Volume name = FDD";

  auto& baseParam = FDDBaseParam::Instance();

  //Rotations used
  //TGeoRotation* Rx90m = new TGeoRotation("Rx90m", 0., -90., 0.);
  //TGeoRotation* Rx90 = new TGeoRotation("Rx90", 0., 90., 0.);
  TGeoRotation* Rx180 = new TGeoRotation("Rx180", 0., 180., 0.);   //   4    |   1
  TGeoRotation* Rz180 = new TGeoRotation("Rz180", 180., 0., 0.);   // --------------->  x
  TGeoRotation* Ry180 = new TGeoRotation("Ry180", 180., 180., 0.); //   3    |   2
  //TGeoRotation* Ry90m = new TGeoRotation("Ry90m", 90., -90., -90.);
  //TGeoRotation* Ry90 = new TGeoRotation("Ry90", 90., 90., -90.);
  //TGeoRotation* Rz90 = new TGeoRotation("Rz90", 90., 0., 0.);

  const Double_t kZbegFrontBar = 1949.1; // Begining of Front Bar

  //Medium for FDA
  TGeoMedium* medFDASci = gGeoManager->GetMedium("FDD_BC420");
  //Medium for FDC
  TGeoMedium* medFDCSci = gGeoManager->GetMedium("FDD_BC420");

  // FDA Scintillator Pad
  const Double_t kFDACellSideY = 21.6;
  const Double_t kFDACellSideX = 18.1;
  // FDC Scintillator Pad
  const Double_t kFDCCellSideY = 21.6;
  const Double_t kFDCCellSideX = 18.1;
  // WLS bar          :  0.40 cm ( 4.0 mm )
  // Wrapping         :  0.20 cm ( 2.0 mm )
  // Aluminnized Mylar:  0.01 cm ( 0.1 mm )
  // Fishing line     :  0.04 cm ( 0.4 mm )
  // total shift on X :  0.65 cm
  // total shift on Y :  0.21 cm
  const Double_t kShiftX = 0.54;
  const Double_t kShiftY = 0.10;
  const Double_t kFDACelldz = 2.54;
  const Double_t kFDCCelldz = 2.54;
  const Double_t kFDABeamPipeR = 6.20; // Radius of beam pipe hole for FDD_A (Diameter  12.4 cm)
  const Double_t kFDCBeamPipeR = 3.70; // Radius of beam pipe hole for FDD_C (Diameter   7.4 cm)
  const Int_t kColorFDA = kGreen;
  const Int_t kColorFDC = kGreen;
  Double_t X = kShiftX + kFDACellSideX * 0.5;
  Double_t Y = kShiftY + kFDACellSideY * 0.5;
  Double_t WLS_dx = 0.4;
  Double_t WLS_dz = 2.5;
  Double_t WLS_SideA_Long_dy = 24.20;  // 24.2;
  Double_t WLS_SideC_Long_dy = 24.20;  // 24.2;
  Double_t WLS_SideA_Short_dy = 18.20; // 18.41;
  Double_t WLS_SideC_Short_dy = 20.70; // 20.91;

  // Creating FDA WLS bars
  TGeoVolume* vFDA_WLS_s = new TGeoVolume("FDAWLSshort", new TGeoBBox("shFDAWLSbarShort", WLS_dx / 2.0, WLS_SideA_Short_dy / 2.0, WLS_dz / 2.0), medFDASci);
  TGeoVolume* vFDA_WLS_l = new TGeoVolume("FDAWLSlong", new TGeoBBox("shFDAWLSbarLong", WLS_dx / 2.0, WLS_SideA_Long_dy / 2.0, WLS_dz / 2.0), medFDASci);
  vFDA_WLS_l->SetLineColor(kRed);
  vFDA_WLS_s->SetLineColor(kRed);

  // Creating FDC WLS bars
  TGeoVolume* vFDC_WLS_s = new TGeoVolume("FDCWLSshort", new TGeoBBox("shFDCWLSbarShort", WLS_dx / 2.0, WLS_SideC_Short_dy / 2.0, WLS_dz / 2.0), medFDCSci);
  TGeoVolume* vFDC_WLS_l = new TGeoVolume("FDCWLSlong", new TGeoBBox("shFDCWLSbarLong", WLS_dx / 2.0, WLS_SideC_Long_dy / 2.0, WLS_dz / 2.0), medFDCSci);
  vFDC_WLS_l->SetLineColor(kRed);
  vFDC_WLS_s->SetLineColor(kRed);

  // Make FDA scintillator pad
  new TGeoBBox("shFDAbox", kFDACellSideX / 2.0, kFDACellSideY / 2.0, kFDACelldz / 2.0);
  new TGeoTube("shFDAHole", 0., kFDABeamPipeR, kFDACelldz);
  (new TGeoTranslation("trFDAbox", X, Y, 0.))->RegisterYourself();

  TGeoVolume* vFDApad = new TGeoVolume("FDApad", new TGeoCompositeShape("shFDApad", "shFDAbox:trFDAbox-shFDAHole"), medFDASci);
  vFDApad->SetLineColor(kColorFDA);

  TGeoVolume* secFDA = new TGeoVolumeAssembly("FDAsec");
  // Add PAD
  Double_t fX_FDA_WLS_s = 0.1 + WLS_dx / 2.0;
  Double_t fX_FDA_WLS_l = kShiftX + WLS_dx / 2.0 + kFDACellSideX + 0.04;
  secFDA->AddNode(vFDApad, 1);
  secFDA->AddNode(vFDA_WLS_s, 1, new TGeoTranslation(fX_FDA_WLS_s, kFDABeamPipeR + WLS_SideA_Short_dy / 2.0, 0.0));
  secFDA->AddNode(vFDA_WLS_l, 1, new TGeoTranslation(fX_FDA_WLS_l, kShiftY + WLS_SideA_Long_dy / 2.0, 0.0));

  /// Assembling FDA adding 4 sectors                                       //  Sectors
  TGeoVolume* vFDAarray = new TGeoVolumeAssembly("FDA"); //        ^ y
  vFDAarray->AddNode(secFDA, 1);                         //        |
  vFDAarray->AddNode(secFDA, 2, Ry180);                  //   2    |   1
  vFDAarray->AddNode(secFDA, 3, Rz180);                  // --------------->  x
  vFDAarray->AddNode(secFDA, 4, Rx180);                  //   3    |   4

  /// FDD_C in the tunnel

  new TGeoBBox("shFDCbox", kFDCCellSideX / 2.0, kFDCCellSideY / 2.0, kFDCCelldz / 2.0);
  new TGeoTube("shFDCHole", 0., kFDCBeamPipeR, kFDCCelldz);
  X = kShiftX + kFDCCellSideX * 0.5;
  Y = kShiftY + kFDCCellSideY * 0.5;
  (new TGeoTranslation("trFDCbox", X, Y, 0.))->RegisterYourself();
  TGeoVolume* vFDCpad = new TGeoVolume("FDCpad", new TGeoCompositeShape("shFDCpad", "shFDCbox:trFDCbox-shFDCHole"), medFDCSci);
  vFDCpad->SetLineColor(kColorFDC);

  /// Creating Sector for Tunnel (Asembly:  Scintillator Pad + Light guide)
  TGeoVolume* voFDC = new TGeoVolumeAssembly("FDCsec");
  // Add PAD
  voFDC->AddNode(vFDCpad, 1);
  // Add FDC WLS Short bar
  voFDC->AddNode(vFDC_WLS_s, 1, new TGeoTranslation(0.1 + WLS_dx / 2.0, kFDCBeamPipeR + WLS_SideC_Short_dy / 2.0, 0.0));
  // Add FDC WLS Long  bar
  voFDC->AddNode(vFDC_WLS_l, 1, new TGeoTranslation(0.04 + WLS_dx / 2.0 + kFDCCellSideX + kShiftX, kShiftY + WLS_SideC_Long_dy / 2.0, 0.0));

  /// Assembling FDC adding the 4 sectors                 //  Sectors
  TGeoVolume* vFDCarray = new TGeoVolumeAssembly("FDC"); //        ^ y
  vFDCarray->AddNode(voFDC, 1);                          //        |
  vFDCarray->AddNode(voFDC, 2, Ry180);                   //   2    |   1
  vFDCarray->AddNode(voFDC, 3, Rz180);                   // --------------->  x
  vFDCarray->AddNode(voFDC, 4, Rx180);                   //   3    |   4
                                                         //        |
  // A side
  TGeoVolume* vCaveRB24 = gGeoManager->GetVolume("caveRB24");
  if (!vCaveRB24) {
    LOG(fatal) << "Could not find the top volume for A-side";
  }
  const Float_t kPosFDA = 1696.67 - 1313.347 - 75.; // z-center of assembly (cm)
  if (baseParam.isExperimental) {
     mChannelCounter = 0;
     float rmin = baseParam.rminA, rmax = baseParam.rmaxA;
     auto vFDDA = buildModuleScint(rmin, rmax);
     vFDDA->SetName("FDDA");
     vCaveRB24->AddNode(vFDDA, 1, new TGeoTranslation(0., 0., Constants::zScintA - Constants::offsetRB24));
     int nsect = Constants::nSectScint;
     for (int i = 0; i < nsect; i ++){
       float rb = 0.5 * (rmin + rmax);
       float theta = (2 * i + 1) * TMath::TwoPi() / nsect / 2;
       auto bundleA = buildFiberBundle();
       vCaveRB24->AddNode(bundleA, 1, new TGeoTranslation(rb * TMath::Cos(theta), rb * TMath::Sin(theta), Constants::zScintA - Constants::offsetRB24 + Constants::dzScint / 2 + Constants::lenFiber / 2 ));
     }
  } else {
    vCaveRB24->AddNode(vFDAarray, 1, new TGeoTranslation(0., 0., kPosFDA - kFDACelldz / 2. - 0.1));
    vCaveRB24->AddNode(vFDAarray, 2, new TGeoTranslation(0., 0., kPosFDA + kFDACelldz / 2. + 0.1));
  }

  // C side
  TGeoVolume* vCave = gGeoManager->GetVolume("cave");
  if (!vCave) {
    LOG(fatal) << "Could not find the top volume for C-side";
  }
  if (baseParam.isExperimental) {
     float rmin = baseParam.rminC, rmax = baseParam.rmaxC;
     auto vFDDC = buildModuleScint(rmin, rmax);
     vFDDC->SetName("FDDC");
     vCave->AddNode(vFDDC, 1, new TGeoTranslation(0., 0., Constants::zScintC));
     int nsect = Constants::nSectScint;
     for (int i = 0; i < nsect; i ++){
       float rb = 0.5 * (rmin + rmax);
       float theta = (2 * i + 1) * TMath::TwoPi() / nsect / 2;
       auto bundleC = buildFiberBundle();
       vCave->AddNode(bundleC, 1, new TGeoTranslation(rb * TMath::Cos(theta), rb * TMath::Sin(theta), Constants::zScintC - Constants::dzScint / 2 - Constants::lenFiber / 2 ));
     }
  } else {
    const Float_t kPosFDC = -kZbegFrontBar - 2. - 3.0 - 0.3; // 3.0 = (5.6 + 0.2 + 0.2)/2.;
    vCave->AddNode(vFDCarray, 1, new TGeoTranslation(0., 0., kPosFDC - kFDCCelldz / 2. - 0.23));
    vCave->AddNode(vFDCarray, 2, new TGeoTranslation(0., 0., kPosFDC + kFDCCelldz / 2. + 0.23));
  }
}


TGeoVolumeAssembly* Geometry::buildModuleScint(float rmin, float rmax)
{
  if (rmin >= rmax) {
    LOG(fatal) << "rmin(" << rmin << ") should be less than rmax(" << rmax << "). Please check your input!";
  }

  auto mod = new TGeoVolumeAssembly("");

  const TGeoMedium* medium = gGeoManager->GetMedium("FDD_Scintillator");

  int nrings = Constants::nRingsScint, nsect = Constants::nSectScint;
  float dz = Constants::dzScint;

  LOG(info) << "Scintillator ring: from " << rmin << " to " << rmax;
  for (int ic = 0; ic < nsect; ic++) {
    int cellId = mChannelCounter++;
    std::string nodeName = "fdd_node" + std::to_string(cellId);
    float phi1 = ic * 360. / nsect;
    float phi2 = (ic + 1) * 360. / nsect;
    auto tbs = new TGeoTubeSeg("tbs", rmin, rmax, dz / 2, phi1, phi2);
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

TGeoVolumeAssembly* Geometry::buildFiberBundle()
{
  auto bundle = new TGeoVolumeAssembly("");
  
  const TGeoMedium* medium = gGeoManager->GetMedium("FDD_Fiber");

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
