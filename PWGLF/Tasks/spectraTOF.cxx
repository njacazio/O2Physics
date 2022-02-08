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
/// \file   spectraTOF.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
///
/// \brief Task for the analysis of the spectra with the TOF detector
///

// O2 includes
#include "ReconstructionDataFormats/Track.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/Core/TrackSelection.h"
#include "Common/Core/TrackSelectionDefaults.h"

using namespace o2;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

// Spectra task
struct tofSpectra {
  static constexpr int Np = 9;
  static constexpr const char* pT[Np] = {"e", "#mu", "#pi", "K", "p", "d", "t", "^{3}He", "#alpha"};
  Configurable<float> cfgNSigmaCut{"cfgNSigmaCut", 3, "Value of the Nsigma cut"};
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8f, "Eta range for tracks"};
  Configurable<int> nBinsP{"nBinsP", 3000, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.01, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 20, "Maximum momentum in range"};
  Configurable<bool> isRun2{"isRun2", false, "Flag to process Run 2 data"};

  // Histograms
  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  static constexpr std::string_view hp[Np] = {"p/El", "p/Mu", "p/Pi",
                                              "p/Ka", "p/Pr", "p/De",
                                              "p/Tr", "p/He", "p/Al"};
  static constexpr std::string_view hpt[Np] = {"pt/El", "pt/Mu", "pt/Pi",
                                               "pt/Ka", "pt/Pr", "pt/De",
                                               "pt/Tr", "pt/He", "pt/Al"};
  static constexpr std::string_view hdcaxy[Np] = {"dcaxy/El", "dcaxy/Mu", "dcaxy/Pi",
                                                  "dcaxy/Ka", "dcaxy/Pr", "dcaxy/De",
                                                  "dcaxy/Tr", "dcaxy/He", "dcaxy/Al"};
  static constexpr std::string_view hdcaz[Np] = {"dcaz/El", "dcaz/Mu", "dcaz/Pi",
                                                 "dcaz/Ka", "dcaz/Pr", "dcaz/De",
                                                 "dcaz/Tr", "dcaz/He", "dcaz/Al"};

  TrackSelection globalTrackswoPrim; // Track without cut for primaries

  void init(o2::framework::InitContext&)
  {
    globalTrackswoPrim = getGlobalTrackSelection();
    globalTrackswoPrim.SetMaxDcaXYPtDep([](float pt) { return 3.f + pt; });
    globalTrackswoPrim.SetRequireGoldenChi2(false);
    if (!isRun2) {
      globalTrackswoPrim.SetTrackType(o2::aod::track::TrackTypeEnum::Track);
    }

    const AxisSpec vtxZAxis{100, -20, 20, "Vtx_{z} (cm)"};
    const AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} (GeV/#it{c})"};
    const AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} (GeV/#it{c})"};

    histos.add("event/vertexz", "", HistType::kTH1F, {vtxZAxis});
    auto h = histos.add<TH1>("evsel", "evsel", HistType::kTH1F, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "posZ passed");
    h = histos.add<TH1>("tracksel", "tracksel", HistType::kTH1F, {{10, 0.5, 10.5}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "Eta passed");
    h->GetXaxis()->SetBinLabel(3, "Quality passed");
    h->GetXaxis()->SetBinLabel(4, "TOF passed");
    histos.add("p/Unselected", "Unselected", kTH1F, {pAxis});
    histos.add("pt/Unselected", "Unselected", kTH1F, {ptAxis});
    for (int i = 0; i < Np; i++) {
      histos.add(hp[i].data(), pT[i], kTH1F, {pAxis});
      histos.add(hpt[i].data(), pT[i], kTH1F, {ptAxis});
    }
    histos.add("electronbeta/hp_El", "", kTH1F, {pAxis});
    histos.add("electronbeta/hpt_El", "", kTH1F, {ptAxis});
    histos.add("electronbeta/hlength_El", ";Track Length (cm);Tracks", kTH1D, {{100, 0, 1000}});
    histos.add("electronbeta/htime_El", ";TOF Time (ns);Tracks", kTH1D, {{1000, 0, 600}});
    histos.add("electronbeta/hp_beta_El", ";#it{p} (GeV/#it{c});#beta - #beta_{e};Tracks", kTH2D, {pAxis, {100, -0.01, 0.01}});
    histos.add("electronbeta/hp_betasigma_El", ";#it{p} (GeV/#it{c});(#beta - #beta_{e})/#sigma;Tracks", kTH2D, {pAxis, {100, -5, 5}});

    // DCAxy
    const AxisSpec dcaXyAxis{600, -3.005, 2.995, "DCA_{xy} (cm)"};
    const AxisSpec dcaZAxis{600, -3.005, 2.995, "DCA_{z} (cm)"};
    for (int i = 0; i < Np; i++) {
      histos.add(hdcaxy[i].data(), pT[i], kTH2F, {ptAxis, dcaXyAxis});
      histos.add(hdcaz[i].data(), pT[i], kTH2F, {ptAxis, dcaZAxis});
    }
  }

  template <PID::ID id, typename T>
  void fillParticleHistos(const T& track)
  {
    const float y = TMath::ASinH(track.pt() / TMath::Sqrt(PID::getMass2(id) + track.pt() * track.pt()) * TMath::SinH(track.eta()));
    if (abs(y) > 0.5) {
      return;
    }
    const auto& nsigma = o2::aod::pidutils::tofNSigma<id>(track);
    if (std::abs(nsigma) < 2) {
      histos.fill(HIST(hdcaxy[id]), track.pt(), track.dcaXY());
      histos.fill(HIST(hdcaz[id]), track.pt(), track.dcaZ());
    }
    if (!track.isGlobalTrack()) {
      return;
    }
    if (abs(nsigma) > cfgNSigmaCut) {
      return;
    }
    histos.fill(HIST(hp[id]), track.p());
    histos.fill(HIST(hpt[id]), track.pt());
  }

  using TrackCandidates = soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksExtended,
                                    aod::pidTOFFullEl, aod::pidTOFFullMu, aod::pidTOFFullPi,
                                    aod::pidTOFFullKa, aod::pidTOFFullPr, aod::pidTOFFullDe,
                                    aod::pidTOFFullTr, aod::pidTOFFullHe, aod::pidTOFFullAl,
                                    aod::pidTOFbeta, aod::TOFSignal,
                                    aod::TrackSelection>;

  void process(aod::Collision const& collision,
               TrackCandidates const& tracks)
  {
    histos.fill(HIST("evsel"), 1);
    if (abs(collision.posZ()) > cfgCutVertex) {
      return;
    }
    histos.fill(HIST("evsel"), 2);
    histos.fill(HIST("event/vertexz"), collision.posZ());

    for (const auto& track : tracks) {
      histos.fill(HIST("tracksel"), 1);
      if (abs(track.eta()) > cfgCutEta) {
        continue;
      }
      histos.fill(HIST("tracksel"), 2);
      if (!globalTrackswoPrim.IsSelected(track)) {
        continue;
      }
      histos.fill(HIST("tracksel"), 3);
      if (!track.hasTOF()) {
        continue;
      }
      histos.fill(HIST("tracksel"), 4);

      histos.fill(HIST("p/Unselected"), track.p());
      histos.fill(HIST("pt/Unselected"), track.pt());

      fillParticleHistos<PID::Electron>(track);
      fillParticleHistos<PID::Muon>(track);
      fillParticleHistos<PID::Pion>(track);
      fillParticleHistos<PID::Kaon>(track);
      fillParticleHistos<PID::Proton>(track);
      fillParticleHistos<PID::Deuteron>(track);
      fillParticleHistos<PID::Triton>(track);
      fillParticleHistos<PID::Helium3>(track);
      fillParticleHistos<PID::Alpha>(track);

      //
      if (TMath::Abs(track.separationbetael() < 1.f)) {
        histos.fill(HIST("electronbeta/hp_El"), track.p());
        histos.fill(HIST("electronbeta/hpt_El"), track.pt());
        histos.fill(HIST("electronbeta/hlength_El"), track.length());
        histos.fill(HIST("electronbeta/htime_El"), track.tofSignal() / 1000);
        histos.fill(HIST("electronbeta/hp_beta_El"), track.p(), track.diffbetael());
        histos.fill(HIST("electronbeta/hp_betasigma_El"), track.p(), track.separationbetael());
      }
    }
  } // end of the process function
};  // end of spectra task

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<tofSpectra>(cfgc)};
}
