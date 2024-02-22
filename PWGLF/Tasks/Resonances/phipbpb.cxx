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
// Phi meson spin alignment task
// sourav.kundu@cern.ch

#include <TH1F.h>
#include <TDirectory.h>
#include <THn.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TFile.h>
#include <TH2F.h>
#include <TLorentzVector.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <Math/Vector4D.h>
#include <cmath>
#include <array>
#include <cstdlib>

#include "PWGLF/DataModel/EPCallibrationTables.h"
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/StepTHn.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Common/DataModel/PIDResponse.h"
#include "Common/DataModel/Multiplicity.h"
#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/trackUtilities.h"
#include "CommonConstants/PhysicsConstants.h"
#include "Common/Core/TrackSelection.h"
#include "Framework/ASoAHelpers.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsParameters/GRPObject.h"
#include "DataFormatsParameters/GRPMagField.h"
#include "CCDB/BasicCCDBManager.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using std::array;
struct phipbpb {

  int mRunNumber;
  int multEstimator;
  float d_bz;
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Service<o2::framework::O2DatabasePDG> pdg;

  // CCDB options
  Configurable<std::string> ccdburl{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<std::string> grpPath{"grpPath", "GLO/GRP/GRP", "Path of the grp file"};
  Configurable<std::string> grpmagPath{"grpmagPath", "GLO/Config/GRPMagField", "CCDB path of the GRPMagField object"};
  Configurable<std::string> lutPath{"lutPath", "GLO/Param/MatLUT", "Path of the Lut parametrization"};
  Configurable<std::string> geoPath{"geoPath", "GLO/Config/GeometryAligned", "Path of the geometry file"};

  // events
  Configurable<float> cfgCutVertex{"cfgCutVertex", 10.0f, "Accepted z-vertex range"};
  // track
  Configurable<float> cfgCutCharge{"cfgCutCharge", 0.0, "cut on Charge"};
  Configurable<float> cfgCutPT{"cfgCutPT", 0.2, "PT cut on daughter track"};
  Configurable<float> cfgCutEta{"cfgCutEta", 0.8, "Eta cut on daughter track"};
  Configurable<float> cfgCutDCAxy{"cfgCutDCAxy", 2.0f, "DCAxy range for tracks"};
  Configurable<float> cfgCutDCAz{"cfgCutDCAz", 2.0f, "DCAz range for tracks"};
  Configurable<float> nsigmaCutTPC{"nsigmacutTPC", 3.0, "Value of the TPC Nsigma cut"};
  Configurable<float> nsigmaCutCombined{"nsigmaCutCombined", 3.0, "Value of the TOF Nsigma cut"};
  Configurable<int> cfgNoMixedEvents{"cfgNoMixedEvents", 1, "Number of mixed events per event"};
  Configurable<int> cfgITScluster{"cfgITScluster", 0, "Number of ITS cluster"};
  Configurable<bool> isDeepAngle{"isDeepAngle", false, "Deep Angle cut"};
  Configurable<double> cfgDeepAngle{"cfgDeepAngle", 0.04, "Deep Angle cut value"};
  Configurable<bool> cfgRemoveOutlier{"Removeoutlier", false, "Additional Event Selection"};
  Configurable<double> cfgEvtSelpar0{"cfgEvtSelpar0", 0.0, "Event selection par0"};
  Configurable<double> cfgEvtSelpar1{"cfgEvtSelpar1", 0.0, "Event selection par1"};
  Configurable<double> cfgEvtSelpar2{"cfgEvtSelpar2", 0.0, "Event selection par2"};
  Configurable<double> cfgEvtSelpar3{"cfgEvtSelpar3", 0.0, "Event selection par3"};

  Filter collisionFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter acceptanceFilter = (nabs(aod::track::eta) < cfgCutEta && nabs(aod::track::pt) > cfgCutPT);
  Filter DCAcutFilter = (nabs(aod::track::dcaXY) < cfgCutDCAxy) && (nabs(aod::track::dcaZ) < cfgCutDCAz);
  Filter PIDcutFilter = nabs(aod::pidtpc::tpcNSigmaKa) < nsigmaCutTPC;

  using EventCandidates = soa::Filtered<soa::Join<aod::Collisions, aod::EvSels, aod::FT0Mults, aod::FV0Mults, aod::TPCMults, aod::CentFV0As, aod::CentFT0Ms, aod::CentFT0Cs, aod::CentFT0As, aod::EPCallibrationTables>>;
  using TrackCandidates = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa>>;

  using EventCandidatesMC = soa::Join<aod::Collisions, aod::EvSels, aod::McCollisionLabels>;
  using TrackCandidatesMC = soa::Filtered<soa::Join<aod::Tracks, aod::TracksExtra, aod::TracksDCA, aod::TrackSelection, aod::pidTPCFullKa, aod::pidTOFFullKa, aod::McTrackLabels>>;

  SliceCache cache;
  Partition<TrackCandidates> posTracks = aod::track::signed1Pt > cfgCutCharge;
  Partition<TrackCandidates> negTracks = aod::track::signed1Pt < cfgCutCharge;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};
  void init(o2::framework::InitContext&)
  {
    histos.add("hFTOMvsTPC", "Mult correlation FT0M vs. TPC", kTH2F, {{600, -0.5f, 59999.5f}, {60, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPC", "Mult correlation FT0C vs. TPC", kTH2F, {{600, -0.5f, 59999.5f}, {60, -0.5f, 5999.5f}});
    histos.add("hFTOAvsTPC", "Mult correlation FT0A vs. TPC", kTH2F, {{600, -0.5f, 59999.5f}, {60, -0.5f, 5999.5f}});
    histos.add("hFTOMvsTPCSelected", "Mult correlation FT0M vs. TPC after selection", kTH2F, {{600, -0.5f, 59999.5f}, {60, -0.5f, 5999.5f}});
    histos.add("hFTOCvsTPCSelected", "Mult correlation FT0C vs. TPC after selection", kTH2F, {{600, -0.5f, 59999.5f}, {60, -0.5f, 5999.5f}});
    histos.add("hFTOAvsTPCSelected", "Mult correlation FT0A vs. TPC after selection", kTH2F, {{600, -0.5f, 59999.5f}, {60, -0.5f, 5999.5f}});
    histos.add("hCentrality", "Centrality distribution", kTH1F, {{201, -0.5, 200.5}});
    histos.add("hVtxZ", "Vertex distribution in Z;Z (cm)", kTH1F, {{400, -20.0, 20.0}});
    histos.add("hEta", "Eta distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaxy", "Dcaxy distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hDcaz", "Dcaz distribution", kTH1F, {{200, -1.0f, 1.0f}});
    histos.add("hNsigmaKaonTPC", "NsigmaKaon TPC distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hNsigmaKaonTOF", "NsigmaKaon TOF distribution", kTH1F, {{200, -10.0f, 10.0f}});
    histos.add("hPsiFT0C", "Psi FT0C", kTH2F, {{111, -0.5, 110.5}, {160, -4.0f, 4.0f}});
    histos.add("h3PhiInvMassSame", "Invariant mass of Phi meson Unlike Sign Same Event", kTH3F, {{111, -0.5, 110.5}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
    histos.add("h3PhiInvMassMixed", "Invariant mass of Phi meson Unlike Sign Mix Event", kTH3F, {{111, -0.5, 110.5}, {200, 0.0f, 20.0f}, {200, 0.9, 1.1}});
  }

  double massKa = o2::constants::physics::MassKPlus;

  template <typename T>
  bool selectionTrack(const T& candidate)
  {
    if (!(candidate.isGlobalTrack() || candidate.isPVContributor() || candidate.itsNCls() > cfgITScluster)) {
      return false;
    }
    return true;
  }
  template <typename T>
  bool selectionPID(const T& candidate)
  {
    if (candidate.hasTOF() && ((candidate.tofNSigmaKa() * candidate.tofNSigmaKa()) + (candidate.tpcNSigmaKa() * candidate.tpcNSigmaKa())) < (nsigmaCutCombined * nsigmaCutCombined)) {
      return true;
    }
    return false;
  }
  // deep angle cut on pair to remove photon conversion
  template <typename T1, typename T2>
  bool selectionPair(const T1& candidate1, const T2& candidate2)
  {
    double pt1, pt2, pz1, pz2, p1, p2, angle;
    pt1 = candidate1.pt();
    pt2 = candidate2.pt();
    pz1 = candidate1.pz();
    pz2 = candidate2.pz();
    p1 = candidate1.p();
    p2 = candidate2.p();
    angle = TMath::ACos((pt1 * pt2 + pz1 * pz2) / (p1 * p2));
    if (isDeepAngle && angle < cfgDeepAngle) {
      return false;
    }
    return true;
  }

  ConfigurableAxis axisVertex{"axisVertex", {20, -10, 10}, "vertex axis for bin"};
  ConfigurableAxis axisMultiplicityClass{"axisMultiplicityClass", {20, 0, 100}, "multiplicity percentile for bin"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {2000, 0, 10000}, "TPC multiplicity  for bin"};

  using BinningTypeVertexContributor = ColumnBinningPolicy<aod::collision::PosZ, aod::cent::CentFT0C>;
  TLorentzVector PhiMother, KaonPlus, KaonMinus;
  void processSameEvent(EventCandidates::iterator const& collision, TrackCandidates const& tracks, aod::BCs const&)
  {
    // auto bc = collision.bc_as<aod::BCsWithTimestamps>(); /// adding timestamp to access magnetic field later
    // initCCDB(bc);

    if (!collision.sel8()) {
      return;
    }
    if (!collision.triggereventep()) {
      return;
    }
    auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision.globalIndex(), cache);
    auto centrality = collision.centFT0C();
    auto multFT0M = collision.multFT0M();
    auto multFT0C = collision.multFT0C();
    auto multFT0A = collision.multFT0A();
    auto multTPC = collision.multTPC();
    histos.fill(HIST("hFTOMvsTPC"), multFT0M, multTPC);
    histos.fill(HIST("hFTOCvsTPC"), multFT0C, multTPC);
    histos.fill(HIST("hFTOCvsTPC"), multFT0A, multTPC);
    if (cfgRemoveOutlier) {
      if ((multTPC > (cfgEvtSelpar0 * multFT0C + cfgEvtSelpar1)) || (multTPC < (cfgEvtSelpar2 * multFT0C + cfgEvtSelpar3))) {
        return;
      }
    }
    histos.fill(HIST("hPsiFT0C"), centrality, collision.psiFT0C());
    histos.fill(HIST("hCentrality"), centrality);
    histos.fill(HIST("hVtxZ"), collision.posZ());
    histos.fill(HIST("hFTOMvsTPCSelected"), multFT0M, multTPC);
    histos.fill(HIST("hFTOCvsTPCSelected"), multFT0C, multTPC);
    histos.fill(HIST("hFTOCvsTPCSelected"), multFT0A, multTPC);
    for (auto track1 : posThisColl) {
      if (!selectionTrack(track1) || !selectionPID(track1)) {
        continue;
      }
      histos.fill(HIST("hEta"), track1.eta());
      histos.fill(HIST("hDcaxy"), track1.dcaXY());
      histos.fill(HIST("hDcaz"), track1.dcaZ());
      histos.fill(HIST("hNsigmaKaonTPC"), track1.tpcNSigmaKa());
      histos.fill(HIST("hNsigmaKaonTOF"), track1.tofNSigmaKa());
      auto track1ID = track1.globalIndex();
      for (auto track2 : negThisColl) {
        if (!selectionTrack(track2) || !selectionPID(track2)) {
          continue;
        }
        auto track2ID = track2.globalIndex();
        if (track2ID == track1ID) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        KaonPlus.SetXYZM(track1.px(), track1.py(), track1.pz(), massKa);
        KaonMinus.SetXYZM(track2.px(), track2.py(), track2.pz(), massKa);
        PhiMother = KaonPlus + KaonMinus;
        histos.fill(HIST("h3PhiInvMassSame"), centrality, PhiMother.Pt(), PhiMother.M());
      }
    }
  }
  PROCESS_SWITCH(phipbpb, processSameEvent, "Process Same event", true);
  void processMixedEvent(EventCandidates const& collisions, TrackCandidates const& tracks)
  {
    BinningTypeVertexContributor binningOnPositions{{axisVertex, axisMultiplicity}, true};
    for (auto const& [collision1, collision2] : o2::soa::selfCombinations(binningOnPositions, cfgNoMixedEvents, -1, collisions, collisions)) {
      if (!collision1.sel8() || !collision2.sel8()) {
        continue;
      }
      auto posThisColl = posTracks->sliceByCached(aod::track::collisionId, collision1.globalIndex(), cache);
      auto negThisColl = negTracks->sliceByCached(aod::track::collisionId, collision2.globalIndex(), cache);

      auto centrality1 = collision1.centFT0C();
      auto multFT0C1 = collision1.multFT0C();
      auto multTPC1 = collision1.multTPC();
      auto multFT0C2 = collision2.multFT0C();
      auto multTPC2 = collision2.multTPC();

      if (cfgRemoveOutlier) {
        if ((multTPC1 > (cfgEvtSelpar0 * multFT0C1 + cfgEvtSelpar1)) || (multTPC1 < (cfgEvtSelpar2 * multFT0C1 + cfgEvtSelpar3))) {
          continue;
        }
        if ((multTPC2 > (cfgEvtSelpar0 * multFT0C2 + cfgEvtSelpar1)) || (multTPC2 < (cfgEvtSelpar2 * multFT0C2 + cfgEvtSelpar3))) {
          continue;
        }
      }

      for (auto& [track1, track2] : o2::soa::combinations(o2::soa::CombinationsFullIndexPolicy(posThisColl, negThisColl))) {
        if (!selectionTrack(track1) || !selectionPID(track1)) {
          continue;
        }
        if (!selectionTrack(track2) || !selectionPID(track2)) {
          continue;
        }
        if (!selectionPair(track1, track2)) {
          continue;
        }
        KaonPlus.SetXYZM(track1.px(), track1.py(), track1.pz(), massKa);
        KaonMinus.SetXYZM(track2.px(), track2.py(), track2.pz(), massKa);
        PhiMother = KaonPlus + KaonMinus;
        histos.fill(HIST("h3PhiInvMassMixed"), centrality1, PhiMother.Pt(), PhiMother.M());
      }
    }
  }
  PROCESS_SWITCH(phipbpb, processMixedEvent, "Process Mixed event", false);
};
WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<phipbpb>(cfgc, TaskName{"phipbpb"})};
}