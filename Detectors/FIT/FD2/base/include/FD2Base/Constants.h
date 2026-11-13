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

/// \file   Constants.h
/// \brief  General constants in FD2

#ifndef ALICEO2_FD2_CONSTANTS_
#define ALICEO2_FD2_CONSTANTS_

namespace o2
{
namespace fd2
{
struct Constants {
  static constexpr unsigned int nSectScint = 8;
  static constexpr unsigned int nRingsScint = 1;

  static constexpr float zScintA = 1700.0f;
  static constexpr float zScintC = -1950.0f;

  static constexpr float offsetRB24 = 1387.5;

  static constexpr float etaMinA = 4.8f;
  static constexpr float etaMaxA = 6.3f;
  static constexpr float etaMinC = -7.0f;
  static constexpr float etaMaxC = -4.9f;

  static constexpr float dzScint = 4.0f;

  static constexpr float rMinFiber = 0.05f;  // fiber radius
  static constexpr float rMaxFiber = 0.1f;  // fiber radius
  static constexpr float lenFiber = 12.0f; // fiber length
};

} // namespace fd2
} // namespace o2
#endif
