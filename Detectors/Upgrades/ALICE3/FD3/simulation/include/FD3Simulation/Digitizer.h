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

#ifndef ALICEO2_FD3_DIGITIZER_H
#define ALICEO2_FD3_DIGITIZER_H

#include "CommonDataFormat/InteractionRecord.h"
#include "DataFormatsFD3/Digit.h"
#include "DataFormatsFD3/ChannelData.h"
#include "DataFormatsFD3/MCLabel.h"
#include "DataFormatsFD3/Triggers.h"
#include "FD3Simulation/Detector.h"
#include "FD3Base/Constants.h"
#include "SimulationDataFormat/MCTruthContainer.h"
#include "FD3Simulation/DigitizationConstants.h"
#include <array>
#include <vector>

namespace o2
{
namespace fd3
{
class Digitizer
{
 private:
  using DP = DigitizationConstants;

 public:
  Digitizer()
    : mTimeStamp(0), mIntRecord(), mEventId(-1), mSrcId(-1), mMCLabels(), mCache(), mPmtChargeVsTime(), mNBins(), mNTimeBinsPerBC(), mPmtResponseGlobalRings(), mPmtResponseTemp(), mLastBCCache(), mCfdStartIndex()
  {
  }

  /// Destructor
  ~Digitizer() = default;

  Digitizer(const Digitizer&) = delete;
  Digitizer& operator=(const Digitizer&) = delete;

  void clear();
  void init();

  void setTimeStamp(long t) { mTimeStamp = t; }
  void setEventId(Int_t id) { mEventId = id; }
  void setSrcId(Int_t id) { mSrcId = id; }
  void setInteractionRecord(const InteractionTimeRecord& ir) { mIntRecord = ir; }

  void process(const std::vector<o2::fd3::Hit>& hits, std::vector<o2::fd3::Digit>& digitsBC,
               std::vector<o2::fd3::ChannelData>& digitsCh, std::vector<o2::fd3::DetTrigInput>& digitsTrig,
               o2::dataformats::MCTruthContainer<o2::fd3::MCLabel>& labels);

  void flush(std::vector<o2::fd3::Digit>& digitsBC,
             std::vector<o2::fd3::ChannelData>& digitsCh,
             std::vector<o2::fd3::DetTrigInput>& digitsTrig,
             o2::dataformats::MCTruthContainer<o2::fd3::MCLabel>& labels);

  const InteractionRecord& getInteractionRecord() const { return mIntRecord; }
  InteractionRecord& getInteractionRecord(InteractionRecord& src) { return mIntRecord; }
  uint32_t getOrbit() const { return mIntRecord.orbit; }
  uint16_t getBC() const { return mIntRecord.bc; }

  using ChannelDigitF = std::vector<float>;

  struct BCCache : public o2::InteractionRecord {
    std::vector<o2::fd3::MCLabel> labels;
    std::array<ChannelDigitF, DP::NCELLSTOT> mPmtChargeVsTime = {};

    void clear()
    {
      for (auto& channel : mPmtChargeVsTime) {
        std::fill(std::begin(channel), std::end(channel), 0.);
      }
      labels.clear();
    }

    BCCache& operator=(const o2::InteractionRecord& ir)
    {
      o2::InteractionRecord::operator=(ir);
      return *this;
    }
    void print() const;
  };

 private:
  static constexpr int BCCacheMin = 0, BCCacheMax = 7, NBC2Cache = 1 + BCCacheMax - BCCacheMin;
  /// Create signal pulse based on MC hit
  /// \param mipFraction Fraction of the MIP energy deposited in the cell
  /// \param parID       Particle ID
  /// \param hitTime     Time of the hit
  /// \param hitR        Length to IP from the position of the hit
  /// \param cachedIR    Cached interaction records
  /// \param nCachedIR   Number of cached interaction records
  /// \param detID       Detector cell ID
  void createPulse(float mipFraction, int parID, const double hitTime, const float hitR,
                   std::array<o2::InteractionRecord, NBC2Cache> const& cachedIR, int nCachedIR, const int detID);

  long mTimeStamp;                  // TF (run) timestamp
  InteractionTimeRecord mIntRecord; // Interaction record (orbit, bc) -> InteractionTimeRecord
  Int_t mEventId;                   // ID of the current event
  Int_t mSrcId;                     // signal, background or QED
  std::deque<fd3::MCLabel> mMCLabels;
  std::deque<BCCache> mCache;

  BCCache& setBCCache(const o2::InteractionRecord& ir);
  BCCache* getBCCache(const o2::InteractionRecord& ir);

  void storeBC(const BCCache& bc,
               std::vector<o2::fd3::Digit>& digitsBC,
               std::vector<o2::fd3::ChannelData>& digitsCh,
               std::vector<o2::fd3::DetTrigInput>& digitsTrig,
               o2::dataformats::MCTruthContainer<o2::fd3::MCLabel>& labels);

  std::array<std::vector<float>, DigitizationConstants::NCELLSTOT> mPmtChargeVsTime; // Charge time series aka analogue signal pulse from PM
  unsigned int mNBins;                                                              //
  unsigned int mNTimeBinsPerBC;
  float mBinSize; // Time width of the pulse bin - HPTDC resolution

  /// vectors to store the PMT signal from cosmic muons
  std::vector<double> mPmtResponseGlobalRings;
  std::vector<double> mPmtResponseTemp;

  /// for CFD
  BCCache mLastBCCache;                                    // buffer for the last BC
  std::array<int, DigitizationConstants::NCELLSTOT> mCfdStartIndex; // start indices for the CFD detector

  /// Internal helper methods related to conversion of energy-deposition into el. signal
  Int_t SimulateLightYield(Int_t pmt, Int_t nPhot) const;
  float SimulateTimeCfd(int& startIndex, const ChannelDigitF& pulseLast, const ChannelDigitF& pulse) const;
  float IntegrateCharge(const ChannelDigitF& pulse) const;

  ClassDefNV(Digitizer, 1);
};
} // namespace fd3
} // namespace o2

#endif
