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
/// \file   qaFakeHits.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2024-04-08
/// \brief  Task to analyze the fraction of the true and fake hits depending on where the fake hits are picked
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

struct QaFakeHits {
  // // Particle information
  // static constexpr int nSpecies = o2::track::PID::NIDs; // One per PDG
  // static constexpr const char* particleTitle[nSpecies] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  // static constexpr int PDGs[nSpecies] = {kElectron, kMuonMinus, kPiPlus, kKPlus, kProton, 1000010020, 1000010030, 1000020030, 1000020040};
  // // Particle only selection
  // Configurable<bool> doEl{"do-el", false, "Flag to run with the PDG code of electrons"};
  // Configurable<bool> doMu{"do-mu", false, "Flag to run with the PDG code of muons"};
  // Configurable<bool> doPi{"do-pi", false, "Flag to run with the PDG code of pions"};
  // Configurable<bool> doKa{"do-ka", false, "Flag to run with the PDG code of kaons"};
  // Configurable<bool> doPr{"do-pr", false, "Flag to run with the PDG code of protons"};
  // Configurable<bool> doDe{"do-de", false, "Flag to run with the PDG code of deuterons"};
  // Configurable<bool> doTr{"do-tr", false, "Flag to run with the PDG code of tritons"};
  // Configurable<bool> doHe{"do-he", false, "Flag to run with the PDG code of helium 3"};
  // Configurable<bool> doAl{"do-al", false, "Flag to run with the PDG code of helium 4"};
  // // Track only selection, options to select only specific tracks
  // Configurable<bool> trackSelection{"trackSelection", true, "Local track selection"};
  // Configurable<int> globalTrackSelection{"globalTrackSelection", 0, "Global track selection: 0 -> No Cut, 1 -> kGlobalTrack, 2 -> kGlobalTrackWoPtEta, 3 -> kGlobalTrackWoDCA, 4 -> kQualityTracks, 5 -> kInAcceptanceTracks, 6 -> custom track cuts via Configurable"};
  // // Event selection
  // Configurable<int> nMinNumberOfContributors{"nMinNumberOfContributors", 2, "Minimum required number of contributors to the primary vertex"};
  // Configurable<float> vertexZMin{"vertex-z-min", -10.f, "Minimum position of the generated vertez in Z (cm)"};
  // Configurable<float> vertexZMax{"vertex-z-max", 10.f, "Maximum position of the generated vertez in Z (cm)"};
  // // Histogram configuration
  // ConfigurableAxis ptBins{"ptBins", {200, 0.f, 5.f}, "Pt binning"};
  // // Histograms
  // HistogramRegistry histosPosPdg{"histosPosPdg", {}, OutputObjHandlingPolicy::AnalysisObject};
  // HistogramRegistry histosNegPdg{"histosNegPdg", {}, OutputObjHandlingPolicy::AnalysisObject};
  // static constexpr int nHistograms = nSpecies * 2;

  // // Pt
  // static constexpr std::string_view hPtIts[nHistograms] = {"MC/el/pos_pdg/pt/its", "MC/mu/pos_pdg/pt/its", "MC/pi/pos_pdg/pt/its",
  //                                                          "MC/ka/pos_pdg/pt/its", "MC/pr/pos_pdg/pt/its", "MC/de/pos_pdg/pt/its",
  //                                                          "MC/tr/pos_pdg/pt/its", "MC/he/pos_pdg/pt/its", "MC/al/pos_pdg/pt/its",
  //                                                          "MC/el/neg_pdg/pt/its", "MC/mu/neg_pdg/pt/its", "MC/pi/neg_pdg/pt/its",
  //                                                          "MC/ka/neg_pdg/pt/its", "MC/pr/neg_pdg/pt/its", "MC/de/neg_pdg/pt/its",
  //                                                          "MC/tr/neg_pdg/pt/its", "MC/he/neg_pdg/pt/its", "MC/al/neg_pdg/pt/its"};

  // static const char* particleName(int pdgSign, o2::track::PID::ID id)
  // {
  //   return Form("%s %s", pdgSign == 0 ? "Positive PDG" : "Negative PDG", o2::track::PID::getName(id));
  // }

  // template <int pdgSign, o2::track::PID::ID id>
  // void makeMCHistograms(const bool doMakeHistograms)
  // {
  //   if (!doMakeHistograms) {
  //     return;
  //   }

  //   if constexpr (pdgSign == 0) {
  //     if (!doPositivePDG) { // Positive
  //       return;
  //     }
  //   } else if constexpr (pdgSign == 1) {
  //     if (!doNegativePDG) { // Negative
  //       return;
  //     }
  //   } else {
  //     LOG(fatal) << "Can't interpret pdgSign " << pdgSign;
  //   }

  //   const AxisSpec axisPt{ptBins, "#it{p}_{T} (GeV/#it{c})"};
  //   const AxisSpec axisP{ptBins, "#it{p} (GeV/#it{c})"};
  //   const AxisSpec axisEta{etaBins, "#it{#eta}"};
  //   const AxisSpec axisY{yBins, "#it{y}"};
  //   const AxisSpec axisPhi{phiBins, "#it{#varphi} (rad)"};

  //   const char* partName = particleName(pdgSign, id);
  //   LOG(info) << "Preparing histograms for particle: " << partName << " pdgSign " << pdgSign;

  //   const TString tagPt = Form("%s #it{#eta} [%.2f,%.2f] #it{y} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
  //                              partName,
  //                              etaMin, etaMax,
  //                              yMin, yMax,
  //                              phiMin, phiMax);

  //   const TString tagEta = Form("%s #it{p}_{T} [%.2f,%.2f] #it{y} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
  //                               partName,
  //                               ptMin, ptMax,
  //                               yMin, yMax,
  //                               phiMin, phiMax);

  //   const TString tagY = Form("%s #it{p}_{T} [%.2f,%.2f] #it{#eta} [%.2f,%.2f] #it{#varphi} [%.2f,%.2f]",
  //                             partName,
  //                             ptMin, ptMax,
  //                             etaMin, etaMax,
  //                             phiMin, phiMax);

  //   const TString tagPhi = Form("%s #it{p}_{T} [%.2f,%.2f] #it{#eta} [%.2f,%.2f] #it{y} [%.2f,%.2f]",
  //                               partName,
  //                               ptMin, ptMax,
  //                               etaMin, etaMax,
  //                               yMin, yMax);

  //   const TString tagPtEta = Form("%s #it{#varphi} [%.2f,%.2f] #it{y} [%.2f,%.2f]",
  //                                 partName,
  //                                 phiMin, phiMax,
  //                                 yMin, yMax);

  //   const int histogramIndex = id + pdgSign * nSpecies;
  //   HistogramRegistry* registry = &histosPosPdg;
  //   if (pdgSign == 1) {
  //     registry = &histosNegPdg;
  //   }

  //   registry->add(hPtIts[histogramIndex].data(), "ITS fakes " + tagPt, kTH1D, {axisPt});

  //   LOG(info) << "Done with particle: " << partName;
  // }

  // void init(InitContext&)
  // {
  //   static_for<0, 1>([&](auto pdgSign) {
  //     makeMCHistograms<pdgSign, o2::track::PID::Electron>(doEl);
  //     makeMCHistograms<pdgSign, o2::track::PID::Muon>(doMu);
  //     makeMCHistograms<pdgSign, o2::track::PID::Pion>(doPi);
  //     makeMCHistograms<pdgSign, o2::track::PID::Kaon>(doKa);
  //     makeMCHistograms<pdgSign, o2::track::PID::Proton>(doPr);
  //     makeMCHistograms<pdgSign, o2::track::PID::Deuteron>(doDe);
  //     makeMCHistograms<pdgSign, o2::track::PID::Triton>(doTr);
  //     makeMCHistograms<pdgSign, o2::track::PID::Helium3>(doHe);
  //     makeMCHistograms<pdgSign, o2::track::PID::Alpha>(doAl);
  //   });
  // }

  // template <int pdgSign, o2::track::PID::ID id, typename particleType>
  // bool isPdgSelected(const particleType& mcParticle)
  // {
  //   static_assert(pdgSign == 0 || pdgSign == 1);
  //   static_assert(id > 0 || id < nSpecies);

  //   // Selecting a specific PDG
  //   if constexpr (pdgSign == 0) {
  //     return mcParticle.pdgCode() == PDGs[id];
  //   } else {
  //     return mcParticle.pdgCode() == -PDGs[id];
  //   }
  // }

  // template <int pdgSign, o2::track::PID::ID id, typename trackType>
  // void fillMCTrackHistograms(const trackType& track, const bool doMakeHistograms)
  // {
  //   static_assert(pdgSign == 0 || pdgSign == 1);
  //   if (!doMakeHistograms) {
  //     return;
  //   }

  //   if constexpr (pdgSign == 0) {
  //     if (!doPositivePDG) {
  //       return;
  //     }
  //   } else {
  //     if (!doNegativePDG) {
  //       return;
  //     }
  //   }

  //   HistogramRegistry* h = &histosPosPdg;
  //   if constexpr (pdgSign == 1) {
  //     h = &histosNegPdg;
  //   }

  //   constexpr int histogramIndex = id + pdgSign * nSpecies;
  //   LOG(debug) << "fillMCTrackHistograms for pdgSign '" << pdgSign << "' and id '" << static_cast<int>(id) << "' " << particleName(pdgSign, id) << " with index " << histogramIndex;
  //   const auto mcParticle = track.mcParticle();

  //   if (!isPdgSelected<pdgSign, id>(mcParticle)) { // Selecting PDG code
  //     return;
  //   }

  //   if (passedITS) {
  //     h->fill(HIST(hPtIts[histogramIndex]), mcParticle.pt());
  //   }
  //   if (passedTPC) {
  //     h->fill(HIST(hPtTpc[histogramIndex]), mcParticle.pt());
  //   }
  // }

  // MC process
  // Preslice<o2::aod::Tracks> perCollision = o2::aod::track::collisionId;
  void process(o2::aod::McCollision const& mcCollision,
              //  o2::soa::SmallGroups<o2::soa::Join<o2::aod::Collisions, o2::aod::McCollisionLabels, o2::aod::EvSels>> const& collisions,
              //  o2::soa::Join<TrackCandidates, o2::aod::McTrackLabels> const& tracks,
               o2::aod::McParticles const& mcParticles)
  {
    // const auto groupedTracks = tracks.sliceBy(perCollision, collision.globalIndex());

    // // Track loop
    // for (const auto& track : groupedTracks) {
    //   static_for<0, 1>([&](auto pdgSign) {
    //     fillMCTrackHistograms<pdgSign, o2::track::PID::Electron>(track, doEl);
    //     fillMCTrackHistograms<pdgSign, o2::track::PID::Muon>(track, doMu);
    //     fillMCTrackHistograms<pdgSign, o2::track::PID::Pion>(track, doPi);
    //     fillMCTrackHistograms<pdgSign, o2::track::PID::Kaon>(track, doKa);
    //     fillMCTrackHistograms<pdgSign, o2::track::PID::Proton>(track, doPr);
    //     fillMCTrackHistograms<pdgSign, o2::track::PID::Deuteron>(track, doDe);
    //     fillMCTrackHistograms<pdgSign, o2::track::PID::Triton>(track, doTr);
    //     fillMCTrackHistograms<pdgSign, o2::track::PID::Helium3>(track, doHe);
    //     fillMCTrackHistograms<pdgSign, o2::track::PID::Alpha>(track, doAl);
    //   });
    // }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<QaFakeHits>(cfgc)};
}
