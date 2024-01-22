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
/// \file   qaPileup.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2024-04-08
/// \brief  Task to analyze the pileup in the simulation
///

// O2 includes
#include "Framework/AnalysisTask.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "ReconstructionDataFormats/DCA.h"
#include "ReconstructionDataFormats/Track.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "PWGLF/DataModel/LFParticleIdentification.h"

using namespace o2::framework;

struct qaPileup {
  HistogramRegistry histos{"histograms", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(InitContext&)
  {
    const AxisSpec axisMask{65536, -0.5f, 65535.5f, "mcMask"};
    histos.add("mcMask", "mcMask", kTH1D, {axisMask});
  }

  //   void process(o2::soa::SmallGroups<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels>>::iterator const& mcCollision,
  void process(o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels>::iterator const& mcCollision,
               o2::aod::McParticles const& mcParticles)
  {
    // LOG(info) << mcCollision.mcMask();
    histos.fill(HIST("mcMask"), mcCollision.mcMask());
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<qaPileup>(cfgc)}; }
