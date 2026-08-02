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
/// \file   dittoEventData.cxx
/// \author Nicolò Jacazio, Università del Piemonte Orientale (IT)
/// \since  2026-08-02
/// \brief  Produce minimal full-phase-space events for DittoMC training
///

#include "ALICE3/DataModel/DittoMc.h"

#include <Framework/AnalysisDataModel.h>
#include <Framework/AnalysisTask.h>
#include <Framework/Configurable.h>
#include <Framework/runDataProcessing.h>

#include <cmath>

using namespace o2;
using namespace o2::framework;

struct DittoEventDataProducer {
  Produces<aod::DittoEvents> outputEvents;
  Produces<aod::DittoParticles> outputParticles;

  Configurable<bool> includeBackgroundParticles{"includeBackgroundParticles",
                                                true,
                                                "Include final-state generator particles tagged as embedded background"};

  Configurable<bool> rejectNonFinite{"rejectNonFinite",
                                     true,
                                     "Reject particles with non-finite four-momentum components"};

  template <typename Particle>
  bool isSelected(Particle const& particle) const
  {
    // DittoMC learns complete generated final states over full phase space.
    if (!particle.producedByGenerator()) {
      return false;
    }
    if (particle.getHepMCStatusCode() != 1) {
      return false;
    }
    if (!includeBackgroundParticles.value && particle.fromBackgroundEvent()) {
      return false;
    }

    if (rejectNonFinite.value &&
        (!std::isfinite(particle.px()) ||
         !std::isfinite(particle.py()) ||
         !std::isfinite(particle.pz()) ||
         !std::isfinite(particle.e()))) {
      return false;
    }

    return true;
  }

  void process(aod::McCollision const& mcCollision,
               aod::McParticles const& mcParticles)
  {
    // Store every MC collision, including events with an empty final state.
    outputEvents(
      static_cast<uint64_t>(mcCollision.globalIndex()),
      mcCollision.weight());

    const auto dittoEventId = outputEvents.lastIndex();

    for (const auto& particle : mcParticles) {
      if (!isSelected(particle)) {
        continue;
      }

      outputParticles(
        dittoEventId,
        particle.pdgCode(),
        particle.px(),
        particle.py(),
        particle.pz(),
        particle.e());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<DittoEventDataProducer>(cfgc)}; }
