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
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "PWGLF/Utils/inelGt.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;

static const std::vector<std::string> parameterNames{"Enable"};
static constexpr int nParameters = 1;
static const int defaultParameters[o2::pwglf::PIDExtended::NIDsTot][nParameters]{{1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {1}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}};
bool enabledArray[o2::pwglf::PIDExtended::NIDsTot];

// Histograms
struct Estimators {
  static const int FT0A = 0;
  static const int FT0C = 1;
  static const int FT0AC = 2;
  static const int FV0A = 3;
  static const int FDDA = 4;
  static const int FDDC = 5;
  static const int FDDAC = 6;
  static const int ZNA = 7;
  static const int ZNC = 8;
  static const int ZEM1 = 9;
  static const int ZEM2 = 10;
  static const int ZPA = 11;
  static const int ZPC = 12;
  static const int ITS = 13;
  static const int nEstimators = 14;

  static constexpr const char* estimatorNames[nEstimators] = {"FT0A",
                                                              "FT0C",
                                                              "FT0AC",
                                                              "FV0A",
                                                              "FDDA",
                                                              "FDDC",
                                                              "FDDAC",
                                                              "ZNA",
                                                              "ZNC",
                                                              "ZEM1",
                                                              "ZEM2",
                                                              "ZPA",
                                                              "ZPC",
                                                              "ITS"};
};

std::array<std::shared_ptr<TH1>, Estimators::nEstimators> hestimators;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsVsITSIB;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsReco;
std::array<std::array<std::shared_ptr<TH2>, o2::pwglf::PIDExtended::NIDsTot>, Estimators::nEstimators> hpt;
std::array<std::array<std::shared_ptr<TH1>, o2::pwglf::PIDExtended::NIDsTot>, Estimators::nEstimators> hyield;

struct mcParticlePrediction {

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosYield{"HistosYield", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosPt{"HistosPt", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis binsEta{"binsEta", {100, -20, 20}, "Binning of the Eta axis"};
  ConfigurableAxis binsPt{"binsPt", {100, 0, 10}, "Binning of the Pt axis"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {1000, 0, 1000}, "Binning of the Multiplicity axis"};
  ConfigurableAxis binsMultiplicityReco{"binsMultiplicityReco", {1000, 0, 10000}, "Binning of the Multiplicity axis"};
  Configurable<LabeledArray<int>> enabledSpecies{"enabledSpecies",
                                                 {defaultParameters[0], o2::pwglf::PIDExtended::NIDsTot, nParameters, o2::pwglf::PIDExtended::arrayNames(), parameterNames},
                                                 "Bethe Bloch parameters"};
  Configurable<bool> selectInelGt0{"selectInelGt0", true, "Select only inelastic events"};
  Service<o2::framework::O2DatabasePDG> pdgDB;

  void init(o2::framework::InitContext&)
  {
    const AxisSpec axisEta{binsEta, "#eta"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisMultiplicity{binsMultiplicity, "Multiplicity"};
    const AxisSpec axisMultiplicityReco{binsMultiplicityReco, "Multiplicity Reco"};

    histos.add("collisions", "collisions", kTH1D, {{10, 0, 10}});
    histos.add("collisionsReco", "collisionsReco", kTH1D, {{10, 0, 10}});
    histos.add("eta/charged", "eta", kTH1D, {axisEta});
    histos.add("eta/neutral", "eta", kTH1D, {axisEta});
    auto h = histos.add<TH1>("particles", "particles", kTH1D, {{o2::pwglf::PIDExtended::NIDsTot, -0.5, -0.5 + o2::pwglf::PIDExtended::NIDsTot}});
    for (int i = 0; i < o2::pwglf::PIDExtended::NIDsTot; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, o2::pwglf::PIDExtended::getName(i));
    }
    h = histos.add<TH1>("particlesPrim", "particlesPrim", kTH1D, {{o2::pwglf::PIDExtended::NIDsTot, -0.5, -0.5 + o2::pwglf::PIDExtended::NIDsTot}});
    for (int i = 0; i < o2::pwglf::PIDExtended::NIDsTot; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, o2::pwglf::PIDExtended::getName(i));
    }

    for (int i = 0; i < Estimators::nEstimators; i++) {
      hestimators[i] = histos.add<TH1>(Form("multiplicity/%s", Estimators::estimatorNames[i]), Estimators::estimatorNames[i], kTH1D, {binsMultiplicity});
      hestimatorsVsITSIB[i] = histos.add<TH2>(Form("multiplicity/vsITSIB/%s", Estimators::estimatorNames[i]), Estimators::estimatorNames[i], kTH2D, {binsMultiplicity, binsMultiplicity});
      hestimatorsReco[i] = histos.add<TH2>(Form("multiplicity/Reco/%s", Estimators::estimatorNames[i]), Estimators::estimatorNames[i], kTH2D, {binsMultiplicity, axisMultiplicityReco});
    }

    for (int i = 0; i < o2::pwglf::PIDExtended::NIDsTot; i++) {
      if (enabledSpecies->get(o2::pwglf::PIDExtended::getName(i), "Enable") != 1) {
        enabledArray[i] = false;
        continue;
      }
      enabledArray[i] = true;
      for (int j = 0; j < Estimators::nEstimators; j++) {
        hpt[j][i] = histosPt.add<TH2>(Form("pt/%s/%s", Estimators::estimatorNames[j], o2::pwglf::PIDExtended::getName(i)), o2::pwglf::PIDExtended::getName(i), kTH2D, {binsPt, axisMultiplicity});
        hyield[j][i] = histosYield.add<TH1>(Form("yield/%s/%s", Estimators::estimatorNames[j], o2::pwglf::PIDExtended::getName(i)), o2::pwglf::PIDExtended::getName(i), kTH1D, {axisMultiplicity});
      }
    }
  }

  template <float etamin, float etamax>
  float countMultInAcceptance(const aod::McParticles& mcParticles)
  {
    static_assert(etamin < etamax, "etamin must be smaller than etamax");
    float counter = 0;
    for (const auto& particle : mcParticles) {

      // primary
      if (!particle.isPhysicalPrimary())
        continue;
      // has pdg
      TParticlePDG* p = pdgDB->GetParticle(particle.pdgCode());
      if (!p)
        continue;
      // is charged
      if (abs(p->Charge()) == 0)
        continue;
      // in acceptance
      if (particle.eta() > etamin && particle.eta() < etamax) {
        counter++;
      }
    }
    return counter;
  }

  template <float etamin, float etamax, bool isNeutral = false>
  float countEnergyInAcceptance(const aod::McParticles& mcParticles)
  {
    static_assert(etamin < etamax, "etamin must be smaller than etamax");
    float counter = 0.f;
    for (const auto& particle : mcParticles) {

      // primary
      if (!particle.isPhysicalPrimary()) {
        continue;
      }
      // has pdg
      TParticlePDG* p = pdgDB->GetParticle(particle.pdgCode());
      if (!p) {
        continue;
      }
      // is neutral
      if constexpr (isNeutral) {
        if (abs(p->Charge()) > 1e-3)
          continue;
      } else {
        if (abs(p->Charge()) <= 1e-3)
          continue;
      }
      // in acceptance
      if (particle.eta() > etamin && particle.eta() < etamax) {
        counter += particle.e();
      }
    }
    return counter;
  }

  float countFT0A(const aod::McParticles& mcParticles) { return countMultInAcceptance<3.5f, 4.9f>(mcParticles); }
  float countFT0C(const aod::McParticles& mcParticles) { return countMultInAcceptance<-3.3f, -2.1f>(mcParticles); }
  float countFV0A(const aod::McParticles& mcParticles) { return countMultInAcceptance<2.2f, 5.1f>(mcParticles); }
  float countFDDA(const aod::McParticles& mcParticles) { return countMultInAcceptance<4.9f, 6.3f>(mcParticles); }
  float countFDDC(const aod::McParticles& mcParticles) { return countMultInAcceptance<-7.f, -4.9f>(mcParticles); }
  float countZNA(const aod::McParticles& mcParticles) { return countEnergyInAcceptance<8.8f, 100.f>(mcParticles); }
  float countZNC(const aod::McParticles& mcParticles) { return countEnergyInAcceptance<-100.f, -8.8f>(mcParticles); }
  float countITSIB(const aod::McParticles& mcParticles) { return countMultInAcceptance<-2.f, 2.f>(mcParticles); }

  void process(aod::McCollision const& mcCollision,
               aod::McParticles const& mcParticles)
  {
    histos.fill(HIST("collisions"), 0.5);
    if (selectInelGt0.value && !o2::pwglf::isINELgt0mc(mcParticles, pdgDB)) {
      return;
    }

    histos.fill(HIST("collisions"), 1.5);
    float nMult[Estimators::nEstimators];
    nMult[Estimators::FT0A] = countFT0A(mcParticles);
    nMult[Estimators::FT0C] = countFT0C(mcParticles);
    nMult[Estimators::FT0AC] = nMult[Estimators::FT0A] + nMult[Estimators::FT0C];
    nMult[Estimators::FV0A] = countFV0A(mcParticles);
    nMult[Estimators::FDDA] = countFDDA(mcParticles);
    nMult[Estimators::FDDC] = countFDDC(mcParticles);
    nMult[Estimators::FDDAC] = nMult[Estimators::FDDA] + nMult[Estimators::FDDC];
    nMult[Estimators::ZNA] = countZNA(mcParticles);
    nMult[Estimators::ZNC] = countZNC(mcParticles);
    nMult[Estimators::ITS] = countITSIB(mcParticles);

    for (int i = 0; i < Estimators::nEstimators; i++) {
      hestimators[i]->Fill(nMult[i]);
      hestimatorsVsITSIB[i]->Fill(nMult[i], nMult[Estimators::ITS]);
    }

    for (const auto& particle : mcParticles) {
      particle.pdgCode();
      const auto id = o2::pwglf::PIDExtended::pdgToId(particle);
      if (id < 0) {
        continue;
      }
      if (!enabledArray[id]) {
        continue;
      }

      if (!particle.isPhysicalPrimary()) {
        continue;
      }

      TParticlePDG* p = pdgDB->GetParticle(particle.pdgCode());
      if (p) {
        if (abs(p->Charge()) > 1e-3)
          histos.fill(HIST("eta/charged"), particle.eta());
        else {
          histos.fill(HIST("eta/neutral"), particle.eta());
        }
      }

      if (abs(particle.y()) > 0.5)
        continue;

      for (int i = 0; i < Estimators::nEstimators; i++) {
        hpt[i][id]->Fill(particle.pt(), nMult[i]);
        hyield[i][id]->Fill(nMult[i]);
      }
    }
  }

  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;
  SliceCache cache;
  void processReco(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::Mults, aod::EvSels>::iterator const& collision,
                   aod::McParticles const& mcParticles,
                   aod::McCollisions const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    if (!collision.sel8()) {
      return;
    }
    const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, collision.mcCollision().globalIndex(), cache);

    histos.fill(HIST("collisionsReco"), 0.5);
    if (selectInelGt0.value && !o2::pwglf::isINELgt0mc(particlesInCollision, pdgDB)) {
      return;
    }

    histos.fill(HIST("collisionsReco"), 1.5);
    float nMult[Estimators::nEstimators];
    nMult[Estimators::FT0A] = countFT0A(particlesInCollision);
    nMult[Estimators::FT0C] = countFT0C(particlesInCollision);
    nMult[Estimators::FT0AC] = nMult[Estimators::FT0A] + nMult[Estimators::FT0C];
    nMult[Estimators::FV0A] = countFV0A(particlesInCollision);
    nMult[Estimators::FDDA] = countFDDA(particlesInCollision);
    nMult[Estimators::FDDC] = countFDDC(particlesInCollision);
    nMult[Estimators::FDDAC] = nMult[Estimators::FDDA] + nMult[Estimators::FDDC];
    nMult[Estimators::ZNA] = countZNA(particlesInCollision);
    nMult[Estimators::ZNC] = countZNC(particlesInCollision);
    nMult[Estimators::ITS] = countITSIB(particlesInCollision);

    hestimatorsReco[Estimators::FT0A]->Fill(nMult[Estimators::FT0A], collision.multFT0A());
    hestimatorsReco[Estimators::FT0C]->Fill(nMult[Estimators::FT0C], collision.multFT0C());
    hestimatorsReco[Estimators::FT0AC]->Fill(nMult[Estimators::FT0AC], collision.multFT0M());
    hestimatorsReco[Estimators::FV0A]->Fill(nMult[Estimators::FV0A], collision.multFV0A());
    hestimatorsReco[Estimators::FDDA]->Fill(nMult[Estimators::FDDA], collision.multFDDA());
    hestimatorsReco[Estimators::FDDC]->Fill(nMult[Estimators::FDDC], collision.multFDDC());
    hestimatorsReco[Estimators::FDDAC]->Fill(nMult[Estimators::FDDAC], collision.multFDDM());
    hestimatorsReco[Estimators::ZNA]->Fill(nMult[Estimators::ZNA], collision.multZNA());
    hestimatorsReco[Estimators::ZNC]->Fill(nMult[Estimators::ZNC], collision.multZNC());
    hestimatorsReco[Estimators::ITS]->Fill(nMult[Estimators::ITS], collision.multNTracksPV());
  }
  PROCESS_SWITCH(mcParticlePrediction, processReco, "Process the reco info", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<mcParticlePrediction>(cfgc)}; }
