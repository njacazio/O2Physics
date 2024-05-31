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
/// \file   mcParticlePrediction.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Francesca Ercolessi francesca.ercolessi@cern.ch
/// \brief Task to build the predictions from the models based on the generated particles
///

#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "PWGLF/Utils/mcParticle.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;

static constexpr const char* chargeNames[2] = {"pos", "neg"};
struct mcParticlePrediction {

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis binsPt{"binsPt", {100, -1, 1}, "Binning of the Pt axis"};
  ConfigurableAxis binsEta{"binsEta", {100, -1, 1}, "Binning of the Eta axis"};
  ConfigurableAxis binsPhi{"binsPhi", {100, -1, 1}, "Binning of the Phi axis"};
  std::array<std::array<std::shared_ptr<TH3>, pid_constants::NIDsTot>, 2> hpt;

  void init(o2::framework::InitContext&)
  {
    histos.add("collisions", "collisions", kTH1D, {{10, 0, 10}});

    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisEta{binsEta, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisPhi{binsPhi, "#it{p}_{T} (GeV/#it{c})"};

    for (int i = 0; i < pid_constants::NIDsTot; i++) {
      for (int j = 0; j < 2; j++) {
        const std::string name = Form("%s/%s", chargeNames[j], PID::getName(i));
        hpt[j][i] = histos.add<TH3>((name + "/pt").c_str(), name.c_str(), kTH3D, {binsPt, binsEta, binsPhi});
      }
    }
  }
  void process(aod::McCollision const& mcCollision,
               aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("collisions"), 0.5);
    for (const auto& particle : mcParticles) {
      const auto id = o2::pwglf::pdgToId(particle);
      if (id == PID::NIDsTot) {
        continue;
      }
      const int chargeIdx = particle.pdgCode() > 0 ? 0 : 1;
      hpt[chargeIdx][id]->Fill(particle.pt(), particle.eta(), particle.phi());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<mcParticlePrediction>(cfgc)}; }
