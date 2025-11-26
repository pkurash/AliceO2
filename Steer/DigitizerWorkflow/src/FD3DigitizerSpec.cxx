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

#include "FD3DigitizerSpec.h"
#include "DataFormatsFD3/ChannelData.h"
#include "DataFormatsFD3/Digit.h"
#include "Framework/ControlService.h"
#include "Framework/ConfigParamRegistry.h"
#include "Framework/DataProcessorSpec.h"
#include "Framework/DataRefUtils.h"
#include "Framework/Lifetime.h"
#include "Headers/DataHeader.h"
#include <TStopwatch.h>
#include "Steer/HitProcessingManager.h" // for DigitizationContext
#include <TChain.h>
#include "SimulationDataFormat/MCTruthContainer.h"
#include "Framework/Task.h"
#include "DataFormatsParameters/GRPObject.h"
#include "FD3Simulation/Digitizer.h"
#include "FD3Simulation/DigitizationConstants.h"
#include "DataFormatsFD3/MCLabel.h"
#include "SimulationDataFormat/MCCompLabel.h"
#include "DetectorsBase/BaseDPLDigitizer.h"
#include "DetectorsRaw/HBFUtils.h"
#include <TFile.h>

using namespace o2::framework;
using SubSpecificationType = o2::framework::DataAllocator::SubSpecificationType;

namespace o2
{
namespace fd3
{

class FD3DPLDigitizerTask : public o2::base::BaseDPLDigitizer
{
  using GRP = o2::parameters::GRPObject;

 public:
  FD3DPLDigitizerTask() : o2::base::BaseDPLDigitizer(), mDigitizer(), mSimChains(), mDigitsCh(), mDigitsBC(), mLabels() {}
  ~FD3DPLDigitizerTask() override = default;

  void initDigitizerTask(framework::InitContext& ic) override
  {
    LOG(debug) << "FD3DPLDigitizerTask:init";
    mDigitizer.init();
    mDisableQED = ic.options().get<bool>("disable-qed"); //TODO: QED implementation to be tested
  }

  void run(framework::ProcessingContext& pc)
  {
    if (mFinished) {
      return;
    }
    LOG(debug) << "FD3DPLDigitizerTask:run";

    // read collision context from input
    auto context = pc.inputs().get<o2::steer::DigitizationContext*>("collisioncontext");
    context->initSimChains(o2::detectors::DetID::FD3, mSimChains);
    const bool withQED = context->isQEDProvided() && !mDisableQED; //TODO: QED implementation to be tested

    mDigitizer.setTimeStamp(context->getGRP().getTimeStart());

    auto& irecords = context->getEventRecords(withQED); //TODO: QED implementation to be tested
    auto& eventParts = context->getEventParts(withQED); //TODO: QED implementation to be tested

    // the interaction record marking the timeframe start
    auto firstTF = InteractionTimeRecord(o2::raw::HBFUtils::Instance().getFirstSampledTFIR(), 0);

    // loop over all composite collisions given from context
    // (aka loop over all the interaction records)
    std::vector<o2::fd3::Hit> hits;
    for (int collID = 0; collID < irecords.size(); ++collID) {
      // Note: Very crude filter to neglect collisions coming before
      // the first interaction record of the timeframe. Remove this, once these collisions can be handled
      // within the digitization routine. Collisions before this timeframe might impact digits of this timeframe.
      // See https://its.cern.ch/jira/browse/O2-5395.
      if (irecords[collID] < firstTF) {
        LOG(info) << "Too early: Not digitizing collision " << collID;
        continue;
      }

      mDigitizer.clear();
      const auto& irec = irecords[collID];
      mDigitizer.setInteractionRecord(irec);
      // for each collision, loop over the constituents event and source IDs
      // (background signal merging is basically taking place here)
      for (auto& part : eventParts[collID]) {
        hits.clear();
        context->retrieveHits(mSimChains, "FD3Hit", part.sourceID, part.entryID, &hits);
        LOG(debug) << "[FD3] For collision " << collID << " eventID " << part.entryID << " found " << hits.size() << " hits ";

        // call actual digitization procedure
        mDigitizer.setEventId(part.entryID);
        mDigitizer.setSrcId(part.sourceID);
        mDigitizer.process(hits, mDigitsBC, mDigitsCh, mDigitsTrig, mLabels);
      }
      LOG(debug) << "[FD3] Has " << mDigitsBC.size() << " BC elements,   " << mDigitsCh.size() << " mDigitsCh elements";
    }

    o2::InteractionTimeRecord terminateIR;
    terminateIR.orbit = 0xffffffff; // supply IR in the infinite future to flush all cached BC
    mDigitizer.setInteractionRecord(terminateIR);
    mDigitizer.flush(mDigitsBC, mDigitsCh, mDigitsTrig, mLabels);

    // here we have all digits and we can send them to consumer (aka snapshot it onto output)
    LOG(info) << "FD3: Sending " << mDigitsBC.size() << " digitsBC and " << mDigitsCh.size() << " digitsCh.";

    // send out to next stage
    pc.outputs().snapshot(Output{"FD3", "DIGITSBC", 0}, mDigitsBC);
    pc.outputs().snapshot(Output{"FD3", "DIGITSCH", 0}, mDigitsCh);
    pc.outputs().snapshot(Output{"FD3", "TRIGGERINPUT", 0}, mDigitsTrig);
    if (pc.outputs().isAllowed({"FD3", "DIGITLBL", 0})) {
      pc.outputs().snapshot(Output{"FD3", "DIGITLBL", 0}, mLabels);
    }
    LOG(info) << "FD3: Sending ROMode= " << mROMode << " to GRPUpdater";
    pc.outputs().snapshot(Output{"FD3", "ROMode", 0}, mROMode);

    // we should be only called once; tell DPL that this process is ready to exit
    pc.services().get<ControlService>().readyToQuit(QuitRequest::Me);
    mFinished = true;
  }

 private:
  bool mFinished = false;
  Digitizer mDigitizer;
  std::vector<TChain*> mSimChains;
  std::vector<o2::fd3::ChannelData> mDigitsCh;
  std::vector<o2::fd3::Digit> mDigitsBC;
  std::vector<o2::fd3::DetTrigInput> mDigitsTrig;
  o2::dataformats::MCTruthContainer<o2::fd3::MCLabel> mLabels; // labels which get filled

  // RS: at the moment using hardcoded flag for continuous readout
  o2::parameters::GRPObject::ROMode mROMode = o2::parameters::GRPObject::ROMode(o2::parameters::GRPObject::CONTINUOUS | o2::parameters::GRPObject::TRIGGERING); // readout mode
  bool mDisableQED = false;
};

o2::framework::DataProcessorSpec getFD3DigitizerSpec(int channel, bool mctruth)
{
  // create the full data processor spec using
  //  a name identifier
  //  input description
  //  algorithmic description (here a lambda getting called once to setup the actual processing function)
  //  options that can be used for this processor (here: input file names where to take the hits)
  std::vector<OutputSpec> outputs;
  outputs.emplace_back("FD3", "DIGITSBC", 0, Lifetime::Timeframe);
  outputs.emplace_back("FD3", "DIGITSCH", 0, Lifetime::Timeframe);
  outputs.emplace_back("FD3", "TRIGGERINPUT", 0, Lifetime::Timeframe);
  if (mctruth) {
    outputs.emplace_back("FD3", "DIGITLBL", 0, Lifetime::Timeframe);
  }
  outputs.emplace_back("FD3", "ROMode", 0, Lifetime::Timeframe);
  outputs.emplace_back("FD3", "DIGITSMCTR", 0, Lifetime::Timeframe);

  return DataProcessorSpec{
    "FD3Digitizer",
    Inputs{InputSpec{"collisioncontext", "SIM", "COLLISIONCONTEXT", static_cast<SubSpecificationType>(channel), Lifetime::Timeframe}},

    outputs,

    AlgorithmSpec{adaptFromTask<FD3DPLDigitizerTask>()},
    Options{{"pileup", VariantType::Int, 1, {"whether to run in continuous time mode"}},
            {"disable-qed", o2::framework::VariantType::Bool, false, {"disable QED handling"}}}};
  //Options{{"pileup", VariantType::Int, 1, {"whether to run in continuous time mode"}}}};
}

} // end namespace fd3
} // end namespace o2
