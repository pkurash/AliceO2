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

#ifndef ALICEO2_FDD_FDDBASEPARAM_
#define ALICEO2_FDD_FDDBASEPARAM_

#include "FDDBase/Constants.h"
#include "FDDBase/Geometry.h"
#include "CommonUtils/ConfigurableParamHelper.h"
#include <string>

namespace o2
{
namespace fdd
{

struct FDDBaseParam : public o2::conf::ConfigurableParamHelper<FDDBaseParam> {

  bool isExperimental = false;

  // Default ring sizes are taken from module positions and coverages
  // zmodA = 1700 cm, 4.8 < eta < 6.3
  // zmodC = -1950 cm, -7.0 < eta < -4.9
  
  float rminA = 6.25;
  float rmaxA = 27.9;

  float rminC = 3.59;
  float rmaxC = 29.34;

  O2ParamDef(FDDBaseParam, "FDDBase");
};

} // namespace fdd
} // namespace o2

#endif
