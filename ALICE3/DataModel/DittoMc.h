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
/// \file   DittoMc.h
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026-08-01
/// \brief  Data to store derived info for ditto training
///

#ifndef ALICE3_DATAMODEL_DITTOMC_H_
#define ALICE3_DATAMODEL_DITTOMC_H_

#include <Framework/ASoA.h>

#include <cstdint>

namespace o2::aod
{

namespace ditto::event
{
DECLARE_SOA_COLUMN(SourceCollisionId, sourceCollisionId, uint64_t); //! Source MC-collision global index
DECLARE_SOA_COLUMN(Weight, weight, float);                          //! MC-event weight
} // namespace ditto::event

DECLARE_SOA_TABLE(DittoEvents, "AOD", "DITTOEVENT",
                  o2::soa::Index<>,
                  ditto::event::SourceCollisionId,
                  ditto::event::Weight);

using DittoEvent = DittoEvents::iterator;

namespace ditto::particle
{
DECLARE_SOA_INDEX_COLUMN(DittoEvent, dittoEvent); //! Parent Ditto event

DECLARE_SOA_COLUMN(PdgCode, pdgCode, int32_t);
DECLARE_SOA_COLUMN(Px, px, float);
DECLARE_SOA_COLUMN(Py, py, float);
DECLARE_SOA_COLUMN(Pz, pz, float);
DECLARE_SOA_COLUMN(E, e, float);
} // namespace ditto::particle

DECLARE_SOA_TABLE(DittoParticles, "AOD", "DITTOPARTICLE",
                  o2::soa::Index<>,
                  ditto::particle::DittoEventId,
                  ditto::particle::PdgCode,
                  ditto::particle::Px,
                  ditto::particle::Py,
                  ditto::particle::Pz,
                  ditto::particle::E);

using DittoParticle = DittoParticles::iterator;
} // namespace o2::aod

#endif // ALICE3_DATAMODEL_DITTOMC_H_
