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

#ifndef ALICEO2_FD3_FD3BASEPARAM_
#define ALICEO2_FD3_FD3BASEPARAM_

#include "FD3Base/GeometryTGeo.h"
#include "FD3Base/Constants.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include <string>

namespace o2
{
namespace fd3
{

//enum FD3GeoVersion = {
//   v1 = 0, 
//   v2, 
//   v3
//};

struct FD3BaseParam : public o2::conf::ConfigurableParamHelper<FD3BaseParam> {

//  float zmodAC = 1500.0f;
//  float zmodAC_extra = 400.0f;

  float dzscint = 4.0f;
  float dzcher  = 4.0f;

  bool modules_extra = false; // switch for modules at z = +-4 m
  
  //int geoVersion = FD3GeoVersion::v1;

  float zscint_v1 = 370.0f;
  float zscint_v2 = 1500.0f;
  float zscint_v3 = 370.0f;

  float zcher_v1 = 430.0f;
  float zcher_v2 = 1560.0f;
  float zcher_v3 = 1560.0f;

  O2ParamDef(FD3BaseParam, "FD3Base");
};

} // namespace fd3
} // namespace o2

#endif
