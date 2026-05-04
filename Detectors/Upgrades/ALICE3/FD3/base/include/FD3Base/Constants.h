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

/// \file   Constants.h
/// \brief  General constants in FV0
///
/// \author Maciej Slupecki, University of Jyvaskyla, Finland

#ifndef ALICEO2_FD3_CONSTANTS_
#define ALICEO2_FD3_CONSTANTS_

namespace o2
{
namespace fd3
{
struct Constants {
  static constexpr unsigned int nsect = 8;
  static constexpr unsigned int nringsScint = 5;
  static constexpr unsigned int nringsCher = 1;

  static constexpr  float dzscint = 4.0f;
  static constexpr  float dzcher  = 4.0f;

  static constexpr float etaMin_scintA_v1 = 2.5f; // modules of eta
  static constexpr float etaMax_scintA_v1 = 5.5f;
  static constexpr float etaMin_scintC_v1 = -4.9f;
  static constexpr float etaMax_scintC_v1 = -2.5f;
  static constexpr float etaMin_cherA_v1  = 4.0f;
  static constexpr float etaMax_cherA_v1  = 5.5f;
  static constexpr float etaMin_cherC_v1  = -4.9f;
  static constexpr float etaMax_cherC_v1  = -4.0f;
  static constexpr float zscintA_v1 = 410.0f;
  static constexpr float zcherA_v1  = 430.0f;
  static constexpr float zscintC_v1 = -410.0f;
  static constexpr float zcherC_v1  = -430.0f;

  static constexpr float etaMin_scint_v2 = 4.0f;
  static constexpr float etaMin_cher_v2  = 5.0f;
  static constexpr float etaMax_scint_v2 = 7.0f;
  static constexpr float etaMax_cher_v2  = 7.0f;
  static constexpr float zscint_v2 = 1500.0f;
  static constexpr float zcher_v2  = 1560.0f;

//  static constexpr float etaMin_scintA_v3 = 2.0f;
//  static constexpr float etaMax_scintA_v3 = 5.0f;
//  static constexpr float etaMin_scintC_v3 = -4.9f;
//  static constexpr float etaMax_scintC_v3 = -2.5f;
//  static constexpr float etaMin_cher_v3  = 5.0f;
//  static constexpr float etaMax_cher_v3  = 7.0f;
//   static constexpr float zscint_v3 = 410.0f;
//   static constexpr float zcher_v3  = 1560.0f;
};

} // namespace fd3
} // namespace o2
#endif
