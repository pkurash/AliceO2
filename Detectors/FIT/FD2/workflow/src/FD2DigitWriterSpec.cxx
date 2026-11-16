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

/// @file   FD2DigitWriterSpec.cxx
#include "FD2Workflow/FD2DigitWriterSpec.h"

namespace o2
{
namespace fd2
{

template <typename T>
using BranchDefinition = framework::MakeRootTreeWriterSpec::BranchDefinition<T>;

o2::framework::DataProcessorSpec getFD2DigitWriterSpec(bool mctruth, bool trigInp)
{
  using InputSpec = framework::InputSpec;
  using MakeRootTreeWriterSpec = framework::MakeRootTreeWriterSpec;
  // Spectators for logging
  auto logger = [](std::vector<o2::fv0::Digit> const& vecDigits) {
    LOG(debug) << "FD2DigitWriter pulled " << vecDigits.size() << " digits";
  };
  // the callback to be set as hook for custom action when the writer is closed
  auto finishWriting = [](TFile* outputfile, TTree* outputtree) {
    const auto* brArr = outputtree->GetListOfBranches();
    int64_t nent = 0;
    for (const auto* brc : *brArr) {
      int64_t n = ((const TBranch*)brc)->GetEntries();
      if (nent && (nent != n)) {
        LOG(error) << "Branches have different number of entries";
      }
      nent = n;
    }
    outputtree->SetEntries(nent);
    outputtree->Write();
    outputfile->Close();
  };

  auto labelsdef = BranchDefinition<o2::dataformats::MCTruthContainer<o2::fv0::MCLabel>>{InputSpec{"labelinput", "FD2", "DIGITSMCTR"},
                                                                                         "FD2DIGITSMCTR", mctruth ? 1 : 0};
  if (trigInp) {
    return MakeRootTreeWriterSpec("FD2DigitWriter",
                                  "fd2digits.root",
                                  "o2sim",
                                  MakeRootTreeWriterSpec::CustomClose(finishWriting),
                                  BranchDefinition<std::vector<o2::fv0::Digit>>{InputSpec{"digitBCinput", "FD2", "DIGITSBC"}, "FD2DigitBC", 1,
                                                                                logger},
                                  BranchDefinition<std::vector<o2::fv0::ChannelData>>{InputSpec{"digitChinput", "FD2", "DIGITSCH"}, "FD2DigitCh"},
                                  BranchDefinition<std::vector<o2::fv0::DetTrigInput>>{InputSpec{"digitTrinput", "FD2", "TRIGGERINPUT"}, "TRIGGERINPUT"},
                                  std::move(labelsdef))();
  } else {
    return MakeRootTreeWriterSpec("FD2DigitWriterRaw",
                                  "o2_fd2digits.root",
                                  "o2sim",
                                  MakeRootTreeWriterSpec::CustomClose(finishWriting),
                                  BranchDefinition<std::vector<o2::fv0::Digit>>{InputSpec{"digitBCinput", "FD2", "DIGITSBC"}, "FD2DigitBC", 1,
                                                                                logger},
                                  BranchDefinition<std::vector<o2::fv0::ChannelData>>{InputSpec{"digitChinput", "FD2", "DIGITSCH"}, "FD2DigitCh"},
                                  std::move(labelsdef))();
  }
}

} // end namespace fd2
} // end namespace o2
