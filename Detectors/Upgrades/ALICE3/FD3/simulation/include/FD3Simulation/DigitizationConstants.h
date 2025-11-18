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

#ifndef ALICEO2_FD3_DIGITIZATION_CONSTANTS
#define ALICEO2_FD3_DIGITIZATION_CONSTANTS

#include "FD3Base/FD3BaseParam.h"
#include "FD3Base/Constants.h"

namespace o2
{
namespace fd3
{
struct DigitizationConstants {
  static constexpr int NCELLSA = Constants::nringsA * Constants::nsect;      // number of scintillator cells side A
  static constexpr int NCELLSC = Constants::nringsC * Constants::nsect;      // number of scintillator cells side C
  static constexpr int NCELLSTOT = NCELLSA + NCELLSC;                        // total number of scintillator cells
  static constexpr float INV_CHARGE_PER_ADC = 1. / 0.6e-12;                  // charge conversion
  static constexpr float TIME_PER_TDCCHANNEL = 0.01302;                      // time conversion from TDC channels to ns
  static constexpr float INV_TIME_PER_TDCCHANNEL = 1. / TIME_PER_TDCCHANNEL; // time conversion from ns to TDC channels
  static constexpr float N_PHOTONS_PER_MEV = 10400;                          // average #photons generated per 1 MeV of deposited energy
};
} // namespace o2
} // namespace fd3
#endif
