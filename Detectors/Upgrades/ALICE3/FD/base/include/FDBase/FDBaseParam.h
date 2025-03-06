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

#ifndef ALICEO2_FD_FDBASEPARAM_
#define ALICEO2_FD_FDBASEPARAM_

#include "FDBase/GeometryTGeo.h"
#include "FDBase/Constants.h"
#include "CommonUtils/ConfigurableParamHelper.h"

namespace o2
{
namespace fd
{
struct FDBaseParam : public o2::conf::ConfigurableParamHelper<FDBaseParam> {

  float zmodA = 1700.0f;
  float zmodC = -1850.0f;
  float dzscint = 4.0f;

  bool withMG = false; // modified geometry

  O2ParamDef(FDBaseParam, "FDBase");
};

} // namespace fd
} // namespace o2

#endif
