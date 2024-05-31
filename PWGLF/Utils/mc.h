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

///
/// \file  mc.h
/// \since 31/05/2024
/// \brief Utilities to handle the MC information
///

#ifndef PWGLF_UTILS_MC_H_
#define PWGLF_UTILS_MC_H_

#include "ReconstructionDataFormats/PID.h"

namespace o2
{
namespace pwglf
{

/// \brief Convert PDG code to PID index
template <typename TrackType>
pid_constants::ID pdgToId(const TrackType& particle)
{
  switch (abs(particle.pdgCode())) {
    case 11:
      return PID::Electron;
    case 13:
      return PID::Muon;
    case 211:
      return PID::Pion;
    case 321:
      return PID::Kaon;
    case 2212:
      return PID::Proton;
    case 1000010020:
      return PID::Deuteron;
    case 1000010030:
      return PID::Triton;
    case 1000020030:
      return PID::Helium3;
    case 1000020040:
      return PID::Alpha;
    case 111:
      return PID::PI0;
    case 22:
      return PID::Photon;
    case 130:
      return PID::K0;
    case 3122:
      return PID::Lambda;
    case 1010010030:
      return PID::HyperTriton;
    case 1010010040:
      return PID::Hyperhydrog4;
    case 3312:
      return PID::XiMinus;
    case 3334:
      return PID::OmegaMinus;
    default:
      LOG(info) << "Cannot identify particle with PDG code " << particle.pdgCode();
      break;
  }
  return PID::NIDsTot;
}

} // namespace pwglf
} // namespace o2

#endif // PWGLF_UTILS_MC_H_
