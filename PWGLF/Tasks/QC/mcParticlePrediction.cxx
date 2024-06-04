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
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "PWGLF/Utils/inelGt.h"

#include "TPDGCode.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::pwglf;

// Particles
static const std::vector<std::string> parameterNames{"Enable"};
static constexpr int nParameters = 1;
static const int defaultParticles[PIDExtended::NIDsTot][nParameters]{{0}, {0}, {1}, {1}, {1}, {0}, {0}, {0}, {0}, {0}, {0}, {1}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}, {0}};
bool enabledParticlesArray[PIDExtended::NIDsTot];

// Estimators
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
  static std::vector<std::string> arrayNames()
  {
    std::vector<std::string> names;
    for (int i = 0; i < nEstimators; i++) {
      names.push_back(estimatorNames[i]);
    }
    return names;
  }
};
bool enabledEstimatorsArray[Estimators::nEstimators];
static const int defaultEstimators[Estimators::nEstimators][nParameters]{{1},  // FT0A
                                                                         {1},  // FT0C
                                                                         {1},  // FT0AC
                                                                         {1},  // FV0A
                                                                         {1},  // FDDA
                                                                         {1},  // FDDC
                                                                         {1},  // FDDAC
                                                                         {0},  // ZNA
                                                                         {0},  // ZNC
                                                                         {0},  // ZEM1
                                                                         {0},  // ZEM2
                                                                         {0},  // ZPA
                                                                         {0},  // ZPC
                                                                         {1}}; // ITS

// Histograms
std::array<std::shared_ptr<TH1>, Estimators::nEstimators> hestimators;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsVsITS;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvGenVsReco;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvRecoVsITS;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvRecoVsFT0A;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvRecoVsBCId;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hestimatorsRecoEvVsBCId;
std::array<std::shared_ptr<TH2>, Estimators::nEstimators> hvertexPosZ;
std::array<std::array<std::shared_ptr<TH2>, PIDExtended::NIDsTot>, Estimators::nEstimators> hpt;
std::array<std::array<std::shared_ptr<TH1>, PIDExtended::NIDsTot>, Estimators::nEstimators> hyield;

struct mcParticlePrediction {

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosRecoEvs{"HistosRecoEvs", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosYield{"HistosYield", {}, OutputObjHandlingPolicy::AnalysisObject};
  HistogramRegistry histosPt{"HistosPt", {}, OutputObjHandlingPolicy::AnalysisObject};
  ConfigurableAxis binsEta{"binsEta", {100, -20, 20}, "Binning of the Eta axis"};
  ConfigurableAxis binsVxy{"binsVxy", {100, -10, 10}, "Binning of the production vertex (x and y) axis"};
  ConfigurableAxis binsVz{"binsVz", {100, -10, 10}, "Binning of the production vertex (z) axis"};
  ConfigurableAxis binsPt{"binsPt", {100, 0, 10}, "Binning of the Pt axis"};
  ConfigurableAxis binsMultiplicity{"binsMultiplicity", {300, -0.5, 299.5}, "Binning of the Multiplicity axis"};
  ConfigurableAxis binsMultiplicityReco{"binsMultiplicityReco", {1000, -0.5, -0.5 + 10000}, "Binning of the Multiplicity axis"};
  Configurable<LabeledArray<int>> enabledSpecies{"enabledSpecies",
                                                 {defaultParticles[0], PIDExtended::NIDsTot, nParameters, PIDExtended::arrayNames(), parameterNames},
                                                 "Particles enabled"};
  Configurable<LabeledArray<int>> enabledEstimators{"enabledEstimators",
                                                    {defaultEstimators[0], Estimators::nEstimators, nParameters, Estimators::arrayNames(), parameterNames},
                                                    "Estimators enabled"};
  Configurable<bool> selectInelGt0{"selectInelGt0", true, "Select only inelastic events"};
  Configurable<bool> discardMismatchedBCs{"discardMismatchedBCs", false, "Select only collisions with matching BC and MC BC"};
  Service<o2::framework::O2DatabasePDG> pdgDB;

  void init(o2::framework::InitContext&)
  {
    const AxisSpec axisEta{binsEta, "#eta"};
    const AxisSpec axisVx{binsVxy, "Vx"};
    const AxisSpec axisVy{binsVxy, "Vy"};
    const AxisSpec axisVz{binsVz, "Vz"};
    const AxisSpec axisPt{binsPt, "#it{p}_{T} (GeV/#it{c})"};
    const AxisSpec axisMultiplicity{binsMultiplicity, "Multiplicity (undefined)"};
    const AxisSpec axisMultiplicityReco{binsMultiplicityReco, "Multiplicity Reco. (undefined)"};

    auto h = histos.add<TH1>("collisions/generated", "collisions", kTH1D, {{10, -0.5, 9.5}});
    h->GetXaxis()->SetBinLabel(1, "Read");
    h->GetXaxis()->SetBinLabel(2, "INELgt0");
    h->GetXaxis()->SetBinLabel(3, "|Z|<10");
    h = histos.add<TH1>("collisions/reconstructed", "collisions", kTH1D, {{10, -0.5, 9.5}});
    h->GetXaxis()->SetBinLabel(1, "Read");
    h->GetXaxis()->SetBinLabel(2, "INELgt0");
    h->GetXaxis()->SetBinLabel(3, "|Z|<10");
    histos.add("particles/eta/charged", "eta", kTH1D, {axisEta});
    histos.add("particles/eta/neutral", "eta", kTH1D, {axisEta});
    histos.add("particles/vtx/x", "Vx", kTH1D, {axisVx});
    histos.add("particles/vtx/y", "Vy", kTH1D, {axisVy});
    histos.add("particles/vtx/z", "Vz", kTH1D, {axisVz});

    for (int i = 0; i < Estimators::nEstimators; i++) {
      if (enabledEstimators->get(Estimators::estimatorNames[i], "Enable") != 1) {
        enabledEstimatorsArray[i] = false;
        continue;
      }
      LOG(info) << "Enabling estimator " << i << " " << Estimators::estimatorNames[i];
      enabledEstimatorsArray[i] = true;
    }

    h = histos.add<TH1>("particles/yields", "particles", kTH1D, {{PIDExtended::NIDsTot, -0.5, -0.5 + PIDExtended::NIDsTot}});
    for (int i = 0; i < PIDExtended::NIDsTot; i++) {
      h->GetXaxis()->SetBinLabel(i + 1, PIDExtended::getName(i));
    }

    for (int i = 0; i < Estimators::nEstimators; i++) {
      if (!enabledEstimatorsArray[i]) {
        continue;
      }
      hestimators[i] = histos.add<TH1>(Form("multiplicity/%s", Estimators::estimatorNames[i]), Estimators::estimatorNames[i], kTH1D, {binsMultiplicity});
      hestimators[i]->GetXaxis()->SetTitle(Form("Multiplicity %s", Estimators::estimatorNames[i]));

      hestimatorsVsITS[i] = histos.add<TH2>(Form("multiplicity/vsITS/%s", Estimators::estimatorNames[i]), Estimators::estimatorNames[i], kTH2D, {binsMultiplicity, binsMultiplicity});
      hestimatorsVsITS[i]->GetXaxis()->SetTitle(Form("Multiplicity %s", Estimators::estimatorNames[i]));
      hestimatorsVsITS[i]->GetYaxis()->SetTitle(Form("Multiplicity %s", Estimators::estimatorNames[Estimators::ITS]));

      hvertexPosZ[i] = histos.add<TH2>(Form("multiplicity/posZ/%s", Estimators::estimatorNames[i]), Estimators::estimatorNames[i], kTH2D, {{200, -20, 20, "pos Z"}, axisMultiplicity});
      hvertexPosZ[i]->GetYaxis()->SetTitle(Form("Multiplicity %s", Estimators::estimatorNames[i]));

      // Reco events
      hestimatorsRecoEvGenVsReco[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/GenVsReco/%s", Estimators::estimatorNames[i]), Estimators::estimatorNames[i], kTH2D, {binsMultiplicity, axisMultiplicityReco});
      hestimatorsRecoEvGenVsReco[i]->GetXaxis()->SetTitle(Form("Multiplicity %s", Estimators::estimatorNames[i]));
      hestimatorsRecoEvGenVsReco[i]->GetYaxis()->SetTitle(Form("Multiplicity Reco. %s", Estimators::estimatorNames[i]));

      hestimatorsRecoEvRecoVsITS[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/RecoVsITS/%s", Estimators::estimatorNames[i]), Estimators::estimatorNames[i], kTH2D, {axisMultiplicityReco, binsMultiplicity});
      hestimatorsRecoEvRecoVsITS[i]->GetXaxis()->SetTitle(Form("Multiplicity Reco. %s", Estimators::estimatorNames[i]));
      hestimatorsRecoEvRecoVsITS[i]->GetYaxis()->SetTitle(Form("Multiplicity %s", Estimators::estimatorNames[Estimators::ITS]));

      hestimatorsRecoEvRecoVsFT0A[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/RecovsFT0A/%s", Estimators::estimatorNames[i]), Estimators::estimatorNames[i], kTH2D, {axisMultiplicityReco, binsMultiplicity});
      hestimatorsRecoEvRecoVsFT0A[i]->GetXaxis()->SetTitle(Form("Multiplicity Reco. %s", Estimators::estimatorNames[i]));
      hestimatorsRecoEvRecoVsFT0A[i]->GetYaxis()->SetTitle(Form("Multiplicity %s", Estimators::estimatorNames[Estimators::FT0A]));

      hestimatorsRecoEvRecoVsBCId[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/RecoVsBCId/%s", Estimators::estimatorNames[i]), Estimators::estimatorNames[i], kTH2D, {{o2::constants::lhc::LHCMaxBunches, -0.5, -0.5 + o2::constants::lhc::LHCMaxBunches}, binsMultiplicityReco});
      hestimatorsRecoEvRecoVsBCId[i]->GetXaxis()->SetTitle(Form("BC in orbit"));
      hestimatorsRecoEvRecoVsBCId[i]->GetYaxis()->SetTitle(Form("Multiplicity Reco. %s", Estimators::estimatorNames[i]));

      hestimatorsRecoEvVsBCId[i] = histosRecoEvs.add<TH2>(Form("multiplicity/Reco/VsBCId/%s", Estimators::estimatorNames[i]), Estimators::estimatorNames[i], kTH2D, {{o2::constants::lhc::LHCMaxBunches, -0.5, -0.5 + o2::constants::lhc::LHCMaxBunches}, binsMultiplicity});
      hestimatorsRecoEvVsBCId[i]->GetXaxis()->SetTitle(Form("BC in orbit"));
      hestimatorsRecoEvVsBCId[i]->GetYaxis()->SetTitle(Form("Multiplicity %s", Estimators::estimatorNames[i]));
    }

    for (int i = 0; i < PIDExtended::NIDsTot; i++) {
      if (enabledSpecies->get(PIDExtended::getName(i), "Enable") != 1) {
        enabledParticlesArray[i] = false;
        continue;
      }
      LOG(info) << "Enabling particle " << i << " " << PIDExtended::getName(i);
      enabledParticlesArray[i] = true;
      for (int j = 0; j < Estimators::nEstimators; j++) {
        if (!enabledEstimatorsArray[j]) {
          continue;
        }
        hpt[j][i] = histosPt.add<TH2>(Form("prediction/pt/%s/%s", Estimators::estimatorNames[j], PIDExtended::getName(i)), PIDExtended::getName(i), kTH2D, {axisPt, axisMultiplicity});
        hpt[j][i]->GetYaxis()->SetTitle(Form("Multiplicity %s", Estimators::estimatorNames[j]));

        hyield[j][i] = histosYield.add<TH1>(Form("prediction/yield/%s/%s", Estimators::estimatorNames[j], PIDExtended::getName(i)), PIDExtended::getName(i), kTH1D, {axisMultiplicity});
        hyield[j][i]->GetYaxis()->SetTitle(Form("Multiplicity %s", Estimators::estimatorNames[j]));
      }
    }
    histos.print();
    histosPt.print();
    histosYield.print();
  }

  template <float etamin, float etamax>
  float countMultInAcceptance(const aod::McParticles& mcParticles)
  {
    static_assert(etamin < etamax, "etamin must be smaller than etamax");
    float counter = 0;
    for (const auto& particle : mcParticles) {

      // primary
      // if (!particle.isPhysicalPrimary()) {
      //   continue;
      // }
      // if ((particle.vx() * particle.vx() + particle.vy() * particle.vy()) > maxProdRadius * maxProdRadius) {
      //   return false;
      // }

      // has pdg
      TParticlePDG* p = pdgDB->GetParticle(particle.pdgCode());
      if (!p) {
        continue;
      }
      // is charged
      if (abs(p->Charge()) == 0) {
        continue;
      }
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
    histos.fill(HIST("collisions/generated"), 0);
    if (selectInelGt0.value && !o2::pwglf::isINELgt0mc(mcParticles, pdgDB)) {
      return;
    }

    histos.fill(HIST("collisions/generated"), 1);
    if (abs(mcCollision.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("collisions/generated"), 2);
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
      if (!enabledEstimatorsArray[i]) {
        continue;
      }

      hestimators[i]->Fill(nMult[i]);
      hestimatorsVsITS[i]->Fill(nMult[i], nMult[Estimators::ITS]);
      hvertexPosZ[i]->Fill(mcCollision.posZ(), nMult[i]);
    }

    for (const auto& particle : mcParticles) {
      particle.pdgCode();
      const auto id = PIDExtended::pdgToId(particle);
      if (id < 0) {
        continue;
      }
      if (!enabledParticlesArray[id]) {
        continue;
      }

      if (!particle.isPhysicalPrimary()) {
        continue;
      }

      TParticlePDG* p = pdgDB->GetParticle(particle.pdgCode());
      if (p) {
        if (abs(p->Charge()) > 1e-3)
          histos.fill(HIST("particles/eta/charged"), particle.eta());
        else {
          histos.fill(HIST("particles/eta/neutral"), particle.eta());
        }
      }

      if (abs(particle.y()) > 0.5) {
        continue;
      }

      histos.fill(HIST("particles/vtx/x"), particle.vx());
      histos.fill(HIST("particles/vtx/y"), particle.vy());
      histos.fill(HIST("particles/vtx/z"), particle.vz() - mcCollision.posZ());

      histos.fill(HIST("particles/yields"), id);
      for (int i = 0; i < Estimators::nEstimators; i++) {
        if (!enabledEstimatorsArray[i]) {
          continue;
        }
        hpt[i][id]->Fill(particle.pt(), nMult[i]);
        hyield[i][id]->Fill(nMult[i]);
      }
    }
  }

  Preslice<aod::McParticles> perMCCol = aod::mcparticle::mcCollisionId;
  SliceCache cache;
  void processReco(soa::Join<aod::Collisions, aod::McCollisionLabels, aod::Mults, aod::EvSels, aod::FT0sCorrected>::iterator const& collision,
                   aod::McParticles const& mcParticles,
                   aod::McCollisions const&,
                   aod::BCs const&)
  {
    if (!collision.has_mcCollision()) {
      return;
    }
    if (!collision.sel8()) {
      return;
    }
    if (!collision.selection_bit(aod::evsel::kIsBBT0A)) {
      return;
    }
    if (!collision.selection_bit(aod::evsel::kIsBBT0C)) {
      return;
    }
    // Check that the BC in data and MC is the same
    if (discardMismatchedBCs.value && collision.bc().globalBC() != collision.mcCollision().bc().globalBC()) {
      return;
    }

    const auto& particlesInCollision = mcParticles.sliceByCached(aod::mcparticle::mcCollisionId, collision.mcCollision().globalIndex(), cache);

    histos.fill(HIST("collisions/reconstructed"), 0);
    if (selectInelGt0.value && !o2::pwglf::isINELgt0mc(particlesInCollision, pdgDB)) {
      return;
    }

    histos.fill(HIST("collisions/reconstructed"), 1);
    if (abs(collision.posZ()) > 10.f) {
      return;
    }
    histos.fill(HIST("collisions/reconstructed"), 2);
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

    float nMultReco[Estimators::nEstimators];
    nMultReco[Estimators::FT0A] = collision.multFT0A();
    nMultReco[Estimators::FT0C] = collision.multFT0C();
    nMultReco[Estimators::FT0AC] = collision.multFT0M();
    nMultReco[Estimators::FV0A] = collision.multFV0A();
    nMultReco[Estimators::FDDA] = collision.multFDDA();
    nMultReco[Estimators::FDDC] = collision.multFDDC();
    nMultReco[Estimators::FDDAC] = collision.multFDDM();
    nMultReco[Estimators::ZNA] = collision.multZNA();
    nMultReco[Estimators::ZNC] = collision.multZNC();
    nMultReco[Estimators::ITS] = collision.multNTracksPV();

    for (int i = 0; i < Estimators::nEstimators; i++) {
      if (!enabledEstimatorsArray[i]) {
        continue;
      }
      hestimatorsRecoEvGenVsReco[i]->Fill(nMult[i], nMultReco[i]);
      hestimatorsRecoEvRecoVsITS[i]->Fill(nMultReco[i], nMult[Estimators::ITS]);
      hestimatorsRecoEvRecoVsFT0A[i]->Fill(nMultReco[i], nMult[Estimators::FT0A]);
      hestimatorsRecoEvRecoVsBCId[i]->Fill(collision.bc().globalBC() % o2::constants::lhc::LHCMaxBunches, nMult[i]);
      hestimatorsRecoEvVsBCId[i]->Fill(collision.bc().globalBC() % o2::constants::lhc::LHCMaxBunches, nMultReco[i]);
    }
  }
  PROCESS_SWITCH(mcParticlePrediction, processReco, "Process the reco info", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<mcParticlePrediction>(cfgc)}; }
