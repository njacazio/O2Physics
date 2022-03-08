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
/// \file   pidTOFBase.h
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \brief  Base to build tasks for TOF PID tasks.
///

#include "Framework/AnalysisTask.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Common/Core/PID/PIDTOF.h"
#include "Common/Core/PID/PIDResponse.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include <CCDB/BasicCCDBManager.h>
#include "TableHelper.h"
#include "TOFBase/EventTimeMaker.h"

using namespace o2;
using namespace o2::pid;
using namespace o2::track;
using namespace o2::framework;
using namespace o2::framework::expressions;

void customize(std::vector<o2::framework::ConfigParamSpec>& workflowOptions)
{
  std::vector<ConfigParamSpec> options{{"add-qa-ev-time", VariantType::Int, 0, {"Produce TOF PID QA histograms for TOF event time"}}};
  std::swap(workflowOptions, options);
}

#include "Framework/runDataProcessing.h"

/// Task to produce the TOF signal from the trackTime information
struct tofSignal {
  o2::framework::Produces<o2::aod::TOFSignal> table;
  bool enableTable = false;

  void init(o2::framework::InitContext& initContext)
  {
    // Checking that the table is requested in the workflow and enabling it
    enableTable = isTableRequiredInWorkflow(initContext, "TOFSignal");
  }
  using Trks = o2::soa::Join<aod::Tracks, aod::TracksExtra>;
  void process(Trks const& tracks)
  {
    if (!enableTable) {
      return;
    }
    table.reserve(tracks.size());
    for (auto& t : tracks) {
      table(o2::pid::tof::TOFSignal<Trks::iterator>::GetTOFSignal(t));
    }
  }
};

/// Selection criteria for tracks used for TOF event time
template <typename trackType>
bool filterForTOFEventTime(const trackType& tr)
{
  return (tr.hasTOF() && tr.p() > 0.5f && tr.p() < 2.f && tr.trackType() == o2::aod::track::TrackTypeEnum::Track);
} // accept all

/// Specialization of TOF event time maker
template <typename trackType,
          bool (*trackFilter)(const trackType&),
          template <typename T, o2::track::PID::ID> typename response,
          typename trackTypeContainer,
          typename responseParametersType>
o2::tof::eventTimeContainer evTimeMakerForTracks(const trackTypeContainer& tracks,
                                                 const responseParametersType& responseParameters,
                                                 const float& diamond = 6.0)
{
  return o2::tof::evTimeMakerFromParam<trackTypeContainer, trackType, trackFilter, response, responseParametersType>(tracks, responseParameters, diamond);
}

/// Task to produce the TOF event time table
struct tofEventTime {
  // Tables to produce
  Produces<o2::aod::TOFEvTime> tableEvTime;
  // Detector response and input parameters
  DetectorResponse response;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> paramfile{"param-file", "", "Path to the parametrization object, if emtpy the parametrization is not taken from file"};
  Configurable<std::string> sigmaname{"param-sigma", "TOFReso", "Name of the parametrization for the expected sigma, used in both file and CCDB mode"};
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> ccdbPath{"ccdbPath", "Analysis/PID/TOF", "Path of the TOF parametrization on the CCDB"};
  Configurable<long> timestamp{"ccdb-timestamp", -1, "timestamp of the object"};

  void init(o2::framework::InitContext& initContext)
  {
    // Getting the parametrization parameters
    ccdb->setURL(url.value);
    ccdb->setTimestamp(timestamp.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    // Not later than now objects
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    //
    const std::vector<float> p = {0.008, 0.008, 0.002, 40.0};
    response.SetParameters(DetectorResponse::kSigma, p);
    const std::string fname = paramfile.value;
    if (!fname.empty()) { // Loading the parametrization from file
      LOG(info) << "Loading exp. sigma parametrization from file" << fname << ", using param: " << sigmaname.value;
      response.LoadParamFromFile(fname.data(), sigmaname.value, DetectorResponse::kSigma);
    } else { // Loading it from CCDB
      std::string path = ccdbPath.value + "/" + sigmaname.value;
      LOG(info) << "Loading exp. sigma parametrization from CCDB, using path: " << path << " for timestamp " << timestamp.value;
      response.LoadParam(DetectorResponse::kSigma, ccdb->getForTimeStamp<Parametrization>(path, timestamp.value));
    }
  }

  using TrksEvTime = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TrackSelection>;
  template <o2::track::PID::ID pid>
  using ResponseImplementationEvTime = o2::pid::tof::ExpTimes<TrksEvTime::iterator, pid>;
  void process(TrksEvTime const& tracks, aod::Collisions const&)
  {
    tableEvTime.reserve(tracks.size());

    int lastCollisionId = -1;      // Last collision ID analysed
    for (auto const& t : tracks) { // Loop on collisions
      if (!t.has_collision()) {    // Track was not assigned, cannot compute event time
        tableEvTime(0.f, 999.f, -1);
      } else if (t.collisionId() == lastCollisionId) { // Event time from this collision is already in the table
        continue;
      }
      /// Create new table for the tracks in a collision
      lastCollisionId = t.collisionId(); /// Cache last collision ID

      const auto tracksInCollision = tracks.sliceBy(aod::track::collisionId, lastCollisionId);
      // First make table for event time
      constexpr float diamond = 6.0;
      const auto evTime = evTimeMakerForTracks<TrksEvTime::iterator, filterForTOFEventTime, o2::pid::tof::ExpTimes>(tracksInCollision, response, diamond);
      static constexpr bool removebias = true;
      int ngoodtracks = 0;
      float et = evTime.eventTime;
      float erret = evTime.eventTimeError;
      for (auto const& trk : tracksInCollision) { // Loop on Tracks
        if constexpr (removebias) {
          o2::tof::removeBias<TrksEvTime::iterator, filterForTOFEventTime>(trk, evTime, ngoodtracks, et, erret, diamond, 2);
        }
        tableEvTime(et, erret, evTime.eventTimeMultiplicity);
      }
    }
  }
};

/// Task that checks the TOF collision time
struct tofPidCollisionTimeQa {
  Configurable<int> nBinsEvTime{"nBinsEvTime", 1000, "Number of bins for the event time"};
  Configurable<float> minEvTime{"minEvTime", -1000.f, "Minimum in range in event time"};
  Configurable<float> maxEvTime{"maxEvTime", 1000.f, "Maximum in range in event time"};
  Configurable<int> nBinsTofSignal{"nBinsTofSignal", 5000, "Number of bins for the tof signal time"};
  Configurable<float> minTofSignal{"minTofSignal", 0.f, "Minimum in range in tof signal time"};
  Configurable<float> maxTofSignal{"maxTofSignal", 100e3, "Maximum in range in tof signal time"};
  Configurable<int> nBinsEvTimeReso{"nBinsEvTimeReso", 1000, "Number of bins in event time resolution"};
  Configurable<float> rangeEvTimeReso{"rangeEvTimeReso", 1000.f, "Range in event time resolution"};
  Configurable<int> nBinsMultiplicity{"nBinsMultiplicity", 1000, "Number of bins for the multiplicity"};
  Configurable<float> rangeMultiplicity{"rangeMultiplicity", 1000.f, "Range for the multiplicity"};
  Configurable<int> logAxis{"logAxis", 0, "Flag to use a log momentum axis"};
  Configurable<int> nBinsP{"nBinsP", 200, "Number of bins for the momentum"};
  Configurable<float> minP{"minP", 0.1f, "Minimum momentum in range"};
  Configurable<float> maxP{"maxP", 5.f, "Maximum momentum in range"};

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::QAObject};
  void init(o2::framework::InitContext& initContext)
  {
    const AxisSpec evTimeAxis{nBinsEvTime, minEvTime, maxEvTime, "TOF event time (ps)"};
    const AxisSpec multAxis{nBinsEvTime, 0, rangeMultiplicity, "Track multiplicity for TOF event time"};
    const AxisSpec evTimeResoAxis{nBinsEvTimeReso, 0, rangeEvTimeReso, "TOF event time resolution (ps)"};
    const AxisSpec tofSignalAxis{nBinsTofSignal, minTofSignal, maxTofSignal, "TOF signal (ps)"};
    AxisSpec pAxis{nBinsP, minP, maxP, "#it{p} GeV/#it{c}"};
    AxisSpec ptAxis{nBinsP, minP, maxP, "#it{p}_{T} GeV/#it{c}"};
    if (logAxis) {
      pAxis.makeLogaritmic();
      ptAxis.makeLogaritmic();
    }
    const AxisSpec collisionAxis{6000, -0.5f, 6000.f - .5f, "Collision index % 6000"};
    const AxisSpec massAxis{1000, 0, 3, "TOF mass (GeV/#it{c}^{2})"};
    const AxisSpec betaAxis{1000, 0, 1.5, "TOF #beta"};
    const AxisSpec deltaAxis{1000, -10000, 10000, "t-t_{ev}-t_{exp}(#pi) (ps)"};
    const AxisSpec lengthAxis{1000, 0, 600, "Track length (cm)"};

    auto h = histos.add<TH1>("eventSelection", "eventSelection", kTH1F, {{10, 0, 10, "Cut passed"}});
    h->GetXaxis()->SetBinLabel(1, "Events read");
    h->GetXaxis()->SetBinLabel(2, "Event selection");
    h->GetXaxis()->SetBinLabel(3, "#sigma_{Ev. time} < 200 ps");
    h->GetXaxis()->SetBinLabel(4, "#sigma_{Ev. time} > 200 ps");
    h = histos.add<TH1>("trackSelection", "trackSelection", kTH1F, {{10, 0, 10, "Cut passed"}});
    h->GetXaxis()->SetBinLabel(1, "Tracks read");
    h->GetXaxis()->SetBinLabel(2, "Track selection");
    h->GetXaxis()->SetBinLabel(3, "hasITS");
    h->GetXaxis()->SetBinLabel(4, "hasTPC");
    h->GetXaxis()->SetBinLabel(5, "hasTOF");
    histos.add("eventTime", "eventTime", kTH1F, {evTimeAxis});
    histos.add("eventTimeReso", "eventTimeReso", kTH1F, {evTimeResoAxis});
    histos.add("eventTimeMult", "eventTimeMult", kTH1F, {multAxis});
    histos.add("eventTimeVsMult", "eventTimeVsMult", kTH2F, {multAxis, evTimeAxis});
    histos.add("eventTimeResoVsMult", "eventTimeResoVsMult", kTH2F, {multAxis, evTimeResoAxis});
    histos.add<TH1>("collisionTime", "collisionTime", kTH1F, {evTimeResoAxis})->GetXaxis()->SetTitle("Collision time (ps)");
    histos.add<TH1>("collisionTimeRes", "collisionTimeRes", kTH1F, {evTimeResoAxis})->GetXaxis()->SetTitle("Collision time resolution (ps)");

    histos.add("tracks/p", "p", kTH1F, {pAxis});
    histos.add("tracks/pt", "pt", kTH1F, {ptAxis});
    histos.add("tracks/length", "length", kTH1F, {lengthAxis});

    histos.add("withtof/p", "p", kTH1F, {pAxis});
    histos.add("withtof/pt", "pt", kTH1F, {ptAxis});
    histos.add("withtof/length", "length", kTH1F, {lengthAxis});
    histos.add("withtof/tofSignal", "tofSignal", kTH1F, {tofSignalAxis});
    histos.add("withtof/beta", "beta", kTH2F, {pAxis, betaAxis});
    histos.add("withtof/delta", "delta", kTH2F, {pAxis, deltaAxis});
    histos.add("withtof/expP", "expP", kTH2F, {pAxis, pAxis});
    histos.add("withtof/mass", "mass", kTH1F, {massAxis});
    histos.add("withtof/tofSignalPerCollision", "tofSignalPerCollision", kTH2S, {collisionAxis, tofSignalAxis});

    histos.add("goodreso/p", "p", kTH1F, {pAxis});
    histos.add("goodreso/pt", "pt", kTH1F, {ptAxis});
    histos.add("goodreso/ptden", "ptden", kTH1F, {ptAxis});
    histos.add("goodreso/length", "length", kTH1F, {lengthAxis});
    histos.add("goodreso/tofSignal", "tofSignal", kTH1F, {tofSignalAxis});
    histos.add("goodreso/beta", "beta", kTH2F, {pAxis, betaAxis});
    histos.add("goodreso/delta", "delta", kTH2F, {pAxis, deltaAxis});
    histos.add("goodreso/expP", "expP", kTH2F, {pAxis, pAxis});
    histos.add("goodreso/mass", "mass", kTH1F, {massAxis});
    histos.add("goodreso/tofSignalPerCollision", "tofSignalPerCollision", kTH2S, {collisionAxis, tofSignalAxis});

    histos.add("badreso/p", "p", kTH1F, {pAxis});
    histos.add("badreso/pt", "pt", kTH1F, {ptAxis});
    histos.add("badreso/ptden", "ptden", kTH1F, {ptAxis});
    histos.add("badreso/length", "length", kTH1F, {lengthAxis});
    histos.add("badreso/tofSignal", "tofSignal", kTH1F, {tofSignalAxis});
    histos.add("badreso/beta", "beta", kTH2F, {pAxis, betaAxis});
    histos.add("badreso/delta", "delta", kTH2F, {pAxis, deltaAxis});
    histos.add("badreso/expP", "expP", kTH2F, {pAxis, pAxis});
    histos.add("badreso/mass", "mass", kTH1F, {massAxis});
    histos.add("badreso/tofSignalPerCollision", "tofSignalPerCollision", kTH2S, {collisionAxis, tofSignalAxis});

    histos.add("goodforevtime/tofSignal", "tofSignal", kTH1F, {tofSignalAxis});
    histos.add("goodforevtime/p", "p", kTH1F, {pAxis});
    histos.add("goodforevtime/pt", "pt", kTH1F, {ptAxis});
    histos.add("goodforevtime/length", "length", kTH1F, {lengthAxis});
    histos.add("goodforevtime/beta", "beta", kTH2F, {pAxis, betaAxis});
    histos.add("goodforevtime/delta", "delta", kTH2F, {pAxis, deltaAxis});
    histos.add("goodforevtime/expP", "expP", kTH2F, {pAxis, pAxis});
    histos.add("goodforevtime/mass", "mass", kTH1F, {massAxis});
    histos.add("goodforevtime/tofSignalPerCollision", "tofSignalPerCollision", kTH2S, {collisionAxis, tofSignalAxis});

    histos.add("withqualitycuts/p", "p", kTH1F, {pAxis});
    histos.add("withqualitycuts/pt", "pt", kTH1F, {ptAxis});
    histos.add("withqualitycuts/length", "length", kTH1F, {lengthAxis});
    histos.add("withqualitycuts/mass", "mass", kTH1F, {massAxis});
  }

  using Trks = soa::Join<aod::Tracks, aod::TracksExtra, aod::TOFSignal, aod::TOFEvTime, aod::TrackSelection>;
  int ncolls = 0;
  void process(aod::Collision const&, Trks const& tracks)
  {
    histos.fill(HIST("eventSelection"), 0.5f);
    histos.fill(HIST("eventSelection"), 1.5f);
    bool eventSet = false;
    for (auto& t : tracks) {
      histos.fill(HIST("trackSelection"), 0.5f);

      if (!t.isGlobalTrack()) {
        continue;
      }
      histos.fill(HIST("trackSelection"), 1.5f);

      if (!t.hasITS()) {
        continue;
      }
      histos.fill(HIST("trackSelection"), 2.5f);
      if (!t.hasTPC()) {
        continue;
      }
      histos.fill(HIST("trackSelection"), 3.5f);

      histos.fill(HIST("tracks/p"), t.p());
      histos.fill(HIST("tracks/pt"), t.pt());
      histos.fill(HIST("tracks/length"), t.length());

      if (t.tofEvTimeErr() > 199.f) {
        histos.fill(HIST("badreso/ptden"), t.pt());
      } else {
        histos.fill(HIST("goodreso/ptden"), t.pt());
      }

      if (!t.hasTOF()) {
        continue;
      }
      histos.fill(HIST("trackSelection"), 4.5f);

      const float beta = o2::pid::tof::Beta<Trks::iterator>::GetBeta(t, t.tofEvTime());
      const float mass = o2::pid::tof::TOFMass<Trks::iterator>::GetTOFMass(t.p(), beta);
      histos.fill(HIST("withtof/p"), t.p());
      histos.fill(HIST("withtof/pt"), t.pt());
      histos.fill(HIST("withtof/length"), t.length());
      histos.fill(HIST("withtof/tofSignal"), t.tofSignal());
      histos.fill(HIST("withtof/beta"), t.p(), beta);
      histos.fill(HIST("withtof/delta"), t.p(), t.tofSignal() - t.tofEvTime() - o2::pid::tof::ExpTimes<Trks::iterator, PID::Pion>::GetExpectedSignal(t));
      histos.fill(HIST("withtof/expP"), t.p(), t.tofExpMom());
      histos.fill(HIST("withtof/mass"), mass);
      histos.fill(HIST("withtof/tofSignalPerCollision"), ncolls % 6000, t.tofSignal());
      if (t.pt() > 0.3 && beta > 0.3) {
        histos.fill(HIST("withqualitycuts/p"), t.p());
        histos.fill(HIST("withqualitycuts/pt"), t.pt());
        histos.fill(HIST("withqualitycuts/length"), t.length());
        histos.fill(HIST("withqualitycuts/mass"), mass);
      }

      if (t.tofEvTimeErr() > 199.f) {
        histos.fill(HIST("badreso/p"), t.p());
        histos.fill(HIST("badreso/pt"), t.pt());
        histos.fill(HIST("badreso/length"), t.length());
        histos.fill(HIST("badreso/tofSignal"), t.tofSignal());
        histos.fill(HIST("badreso/beta"), t.p(), beta);
        histos.fill(HIST("badreso/delta"), t.p(), t.tofSignal() - t.tofEvTime() - o2::pid::tof::ExpTimes<Trks::iterator, PID::Pion>::GetExpectedSignal(t));
        histos.fill(HIST("badreso/expP"), t.p(), t.tofExpMom());
        histos.fill(HIST("badreso/mass"), mass);
        histos.fill(HIST("badreso/tofSignalPerCollision"), ncolls % 6000, t.tofSignal());
      } else {
        histos.fill(HIST("goodreso/p"), t.p());
        histos.fill(HIST("goodreso/pt"), t.pt());
        histos.fill(HIST("goodreso/length"), t.length());
        histos.fill(HIST("goodreso/tofSignal"), t.tofSignal());
        histos.fill(HIST("goodreso/beta"), t.p(), beta);
        histos.fill(HIST("goodreso/delta"), t.p(), t.tofSignal() - t.tofEvTime() - o2::pid::tof::ExpTimes<Trks::iterator, PID::Pion>::GetExpectedSignal(t));
        histos.fill(HIST("goodreso/expP"), t.p(), t.tofExpMom());
        histos.fill(HIST("goodreso/mass"), mass);
        histos.fill(HIST("goodreso/tofSignalPerCollision"), ncolls % 6000, t.tofSignal());
      }

      if (!eventSet) {
        if (t.tofEvTimeErr() > 199.f) {
          histos.fill(HIST("eventSelection"), 2.5f);
        } else {
          histos.fill(HIST("eventSelection"), 3.5f);
        }
        histos.fill(HIST("eventTime"), t.tofEvTime());
        histos.fill(HIST("eventTimeReso"), t.tofEvTimeErr());
        histos.fill(HIST("eventTimeMult"), t.tofEvTimeMult());
        histos.fill(HIST("eventTimeVsMult"), t.tofEvTimeMult(), t.tofEvTime());
        histos.fill(HIST("eventTimeResoVsMult"), t.tofEvTimeMult(), t.tofEvTimeErr());

        histos.fill(HIST("collisionTime"), t.collision().collisionTime());
        histos.fill(HIST("collisionTimeRes"), t.collision().collisionTimeRes());
        eventSet = true;
        ncolls++;
      }
      if (!filterForTOFEventTime(t)) {
        continue;
      }
      histos.fill(HIST("goodforevtime/p"), t.p());
      histos.fill(HIST("goodforevtime/pt"), t.pt());
      histos.fill(HIST("goodforevtime/length"), t.length());
      histos.fill(HIST("goodforevtime/tofSignal"), t.tofSignal());
      histos.fill(HIST("goodforevtime/beta"), t.p(), beta);
      histos.fill(HIST("goodforevtime/delta"), t.p(), t.tofSignal() - t.tofEvTime() - o2::pid::tof::ExpTimes<Trks::iterator, PID::Pion>::GetExpectedSignal(t));
      histos.fill(HIST("goodforevtime/expP"), t.p(), t.tofExpMom());
      histos.fill(HIST("goodforevtime/mass"), mass);
      histos.fill(HIST("goodforevtime/tofSignalPerCollision"), ncolls % 6000, t.tofSignal());
    }
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  auto workflow = WorkflowSpec{adaptAnalysisTask<tofSignal>(cfgc),
                               adaptAnalysisTask<tofEventTime>(cfgc)};
  if (cfgc.options().get<int>("add-qa-ev-time")) {
    workflow.push_back(adaptAnalysisTask<tofPidCollisionTimeQa>(cfgc));
  }
  return workflow;
}
