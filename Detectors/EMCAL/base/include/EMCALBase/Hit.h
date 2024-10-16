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

#ifndef ALICEO2_EMCAL_HIT_H
#define ALICEO2_EMCAL_HIT_H

#include "SimulationDataFormat/BaseHits.h"
#include "CommonUtils/ShmAllocator.h"

namespace o2
{
namespace emcal
{
/// \class Hit
/// \brief EMCAL simulation hit information
/// \ingroup EMCALbase
/// \author Markus Fasel <markus.fasel@cern.ch>, Oak Ridge National Laboratory
/// \since August 31st, 2017
class Hit : public o2::BasicXYZEHit<float>
{
 public:
  /// \brief Default constructor
  Hit() = default;

  /// \brief Hit constructor
  ///
  /// Fully defining information of the EMCAL point (position,
  /// momentum, energy, track, ...)
  ///
  /// \param primary Number of primary particle
  /// \param trackID Index of the track, defined as parent track entering teh EMCAL
  /// \param detID ID of the detector segment
  /// \param initialEnergy Energy of the primary particle enering the EMCAL
  /// \param pos Position vector of the point
  /// \param mom Momentum vector for the particle at the point
  /// \param tof Time of the hit
  /// \param eLoss Energy loss
  Hit(Int_t primary, Int_t trackID, Int_t detID, Double_t initialEnergy, const math_utils::Point3D<float>& pos,
      const math_utils::Vector3D<float>& mom, Double_t tof, Double_t eLoss)
    : o2::BasicXYZEHit<float>(pos.X(), pos.Y(), pos.Z(), tof, eLoss, trackID, detID),
      mPvector(mom),
      mPrimary(primary),
      mInitialEnergy(initialEnergy)
  {
  }

  /// \brief Check whether the points are from the same parent and in the same detector volume
  /// \return True if points are the same (origin and detector), false otherwise
  Bool_t operator==(const Hit& rhs) const;

  /// \brief Sorting points according to parent particle and detector volume
  /// \return True if this point is smaller, false otherwise
  Bool_t operator<(const Hit& rhs) const;

  /// \brief Adds energy loss from the other point to this point
  /// \param rhs EMCAL point to add to this point
  /// \return This point with the summed energy loss
  Hit& operator+=(const Hit& rhs);

  /// \brief Destructor
  ~Hit() = default;

  /// \brief Get the initial energy of the primary particle entering EMCAL
  /// \return Energy of the primary particle entering EMCAL
  Double_t GetInitialEnergy() const { return mInitialEnergy; }

  /// \brief Get Primary particles at the origin of the hit
  /// \return Primary particles at the origin of the hit
  Int_t GetPrimary() const { return mPrimary; }

  /// \brief Set initial energy of the primary particle entering EMCAL
  /// \param energy Energy of the primary particle entering EMCAL
  void SetInitialEnergy(Double_t energy) { mInitialEnergy = energy; }

  /// \brief Set primary particles at the origin of the hit
  /// \param primary Primary particles at the origin of the hit
  void SetPrimary(Int_t primary) { mPrimary = primary; }

  /// \brief Writing point information to an output stream;
  /// \param stream target output stream
  void PrintStream(std::ostream& stream) const;

 private:
  math_utils::Vector3D<float> mPvector; ///< Momentum Vector
  Int_t mPrimary;                       ///< Primary particles at the origin of the hit
  Double32_t mInitialEnergy;            ///< Energy of the parent particle that entered the EMCAL

  ClassDefNV(Hit, 1);
};

/// \brief Creates a new point base on this point but adding the energy loss of the right hand side
/// \param lhs Left hand side of the sum
/// \param rhs Right hand side of the sum
/// \return New EMAL point base on this point
Hit operator+(const Hit& lhs, const Hit& rhs);

/// \brief Output stream operator for EMCAL hits
/// \param stream Stream to write on
/// \param point Hit to be printed
/// \return Stream after printing
std::ostream& operator<<(std::ostream& stream, const Hit& point);
} // namespace emcal
} // namespace o2

#ifdef USESHM
namespace std
{
template <>
class allocator<o2::emcal::Hit> : public o2::utils::ShmAllocator<o2::emcal::Hit>
{
};
} // namespace std
#endif

#endif /* Point_h */
