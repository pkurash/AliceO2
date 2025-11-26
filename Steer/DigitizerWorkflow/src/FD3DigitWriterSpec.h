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

#ifndef STEER_DIGITIZERWORKFLOW_FD3DIGITWRITER_H_
#define STEER_DIGITIZERWORKFLOW_FD3DIGITWRITER_H_

#include "Framework/DataProcessorSpec.h"
#include "DPLUtils/MakeRootTreeWriterSpec.h"
#include "Framework/InputSpec.h"
#include "DataFormatsFD3/ChannelData.h"
#include "DataFormatsFD3/Digit.h"
#include "DataFormatsFD3/MCLabel.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "SimulationDataFormat/MCCompLabel.h"

namespace o2
{
namespace fd3
{

template <typename T>
using BranchDefinition = framework::MakeRootTreeWriterSpec::BranchDefinition<T>;

o2::framework::DataProcessorSpec getFD3DigitWriterSpec(bool mctruth = true)
{
  using InputSpec = framework::InputSpec;
  using MakeRootTreeWriterSpec = framework::MakeRootTreeWriterSpec;
  return MakeRootTreeWriterSpec("FD3DigitWriter",
                                "fd3digits.root",
                                "o2sim",
                                1,
                                BranchDefinition<std::vector<o2::fd3::Digit>>{InputSpec{"fd3digitBCinput", "FD3", "DIGITSBC"}, "FD3DigitBC"},
                                BranchDefinition<std::vector<o2::fd3::ChannelData>>{InputSpec{"fd3digitChinput", "FD3", "DIGITSCH"}, "FD3DigitCh"},
                                BranchDefinition<std::vector<o2::fd3::DetTrigInput>>{InputSpec{"fd3digitTrinput", "FD3", "TRIGGERINPUT"}, "TRIGGERINPUT"},
                                BranchDefinition<o2::dataformats::MCTruthContainer<o2::ft0::MCLabel>>{InputSpec{"FT0labelinput", "FT0", "DIGITSMCTR"}, "FT0DIGITSMCTR", mctruth ? 1 : 0},
                                BranchDefinition<o2::dataformats::MCTruthContainer<o2::fd3::MCLabel>>{InputSpec{"fd3labelinput", "FD3", "DIGITLBL"}, "FD3DigitLabels", mctruth ? 1 : 0})();
}

} // namespace fd3
} // end namespace o2

#endif /* STEER_DIGITIZERWORKFLOW_FD3DIGITWRITER_H_ */
