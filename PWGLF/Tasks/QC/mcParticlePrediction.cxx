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
/// \author Francesca Ercolessi francesca.ercolessi@cern.ch
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
///
/// \brief Task to build the predictions from the models based on the generated particles
///

#include "Framework/HistogramRegistry.h"
#include "Framework/StaticFor.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "ReconstructionDataFormats/PID.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

static constexpr const char* chargeNames[2] = {"pos", "neg"};
struct mcParticlePrediction {

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis binsPt{"binsPt", {100, -1, 1}, "Binning of the Pt axis"};
  std::array<std::array<std::shared_ptr<TH3>, pid_constants::NIDsTot>, 2> hpt;

  pid_constants::ID pdgToId(const aod::McParticle& particle)
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
        LOG(info) << "Cannot identify particle with PDG code " << particle.pdgCode() << " and mass " << particle.mass() << " GeV/c^2";
        break;
    }
    return PID::NIDsTot;
  }

  void init(o2::framework::InitContext&)
  {
    const AxisSpec ptAxis{binsPt, "#it{p}_{T} (GeV/#it{c})"};

    for (int i = 0; i < pid_constants::NIDsTot; i++) {
      for (int j = 0; j < 2; j++) {
        const std::string name = Form("%s/%s", chargeNames[j], pid_constants::getName(i));
        hpt[j][i] = histos.make<TH3D>(name + "/pt", name, binsPt, binsPt, binsPt);
      }
    }
  }
  void process(aod::McCollision const& mcCollision,
               aod::McParticles const& mcParticles)
  {
    for (const auto& particle : mcParticles) {
      const auto id = pdgToId(particle);
      if (id == PID::NIDsTot) {
        continue;
      }
      const int chargeIdx = particle.pdgCode() > 0 ? 0 : 1;
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<mcParticlePrediction>(cfgc)}; }
