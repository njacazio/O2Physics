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
/// \file   mcCentrality.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \author Francesca Ercolessi francesca.ercolessi@cern.ch
/// \since  2024-06-05
/// \brief  Task to produce the table for the equalized multiplicity into centrality bins
///

// O2 includes
#include "CCDB/BasicCCDBManager.h"
#include "ReconstructionDataFormats/Track.h"
#include "CCDB/CcdbApi.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/StaticFor.h"
#include "TableHelper.h"
#include "Framework/O2DatabasePDGPlugin.h"
#include "Common/DataModel/Centrality.h"

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Task to produce the response table
struct mcCentrality {

  // Tables to produce
  Produces<aod::CentRun2V0Ms> centRun2V0M;
  Produces<aod::CentRun2V0As> centRun2V0A;
  Produces<aod::CentRun2SPDTrks> centRun2SPDTracklets;
  Produces<aod::CentRun2SPDClss> centRun2SPDClusters;
  Produces<aod::CentRun2CL0s> centRun2CL0;
  Produces<aod::CentRun2CL1s> centRun2CL1;
  Produces<aod::CentFV0As> centFV0A;
  Produces<aod::CentFT0Ms> centFT0M;
  Produces<aod::CentFT0As> centFT0A;
  Produces<aod::CentFT0Cs> centFT0C;
  Produces<aod::CentFDDMs> centFDDM;
  Produces<aod::CentNTPVs> centNTPV;

  // Input parameters
  Service<o2::ccdb::BasicCCDBManager> ccdb;
  Configurable<std::string> url{"ccdb-url", "http://alice-ccdb.cern.ch", "url of the ccdb repository"};
  Configurable<int64_t> ccdbTimestamp{"ccdb-timestamp", -1, "timestamp of the object used to query in CCDB the detector response. If 0 the object corresponding to the run number is used, if < 0 the latest object is used"};
  Configurable<std::string> path{"path", "/tmp/InputCalibMC.root", "path to calib file or ccdb path if begins with ccdb://"};
  Configurable<bool> selectPrimaries{"selectPrimaries", true, "Select only primary particles"};
  Service<o2::framework::O2DatabasePDG> pdgDB;

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  TH1F* h1dFT0M;
  /*TH1F* h1dFT0A;
  TH1F* h1dFT0C;
  TH1F* h1dFDD;
  TH1F* h1dNTP;*/

  void init(o2::framework::InitContext& initContext)
  {
    // Set up the CCDB
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    ccdb->setFatalWhenNull(false);

    TList* lOfInput;
    if (path.value.rfind("ccdb://", 0) == 0) { // Getting post calib. from CCDB
      path.value.replace(0, 7, "");
      lOfInput = ccdb->get<TList>(path);
      if (!lOfInput) {
        LOG(fatal) << "Could not find the calibration TList from CCDB in path " << path;
        return;
      }
    } else { // Getting post calib. from file
      TFile* f = TFile::Open(path.value.c_str(), "READ");
      if (!f) {
        LOG(fatal) << "The input file " << path << " is not valid";
      }
      if (!f->IsOpen()) {
        LOG(fatal) << "The input file " << f->GetName() << " is not open";
      }
      lOfInput = static_cast<TList*>(f->Get("ccdb_object"));
      if (!lOfInput) {
        f->ls();
        LOG(fatal) << "The input file " << path.value << " does not contain the TList ccdb_object";
      }
    }
    h1dFT0M = static_cast<TH1F*>(lOfInput->FindObject("h1dFT0M"));
    if (!h1dFT0M) {
      lOfInput->ls();
      LOG(fatal) << "Could not open histogram h1dFT0M from TList";
      return;
    }
  }

  template <float etamin, float etamax>
  float countMultInAcceptance(const aod::McParticles& mcParticles)
  {
    static_assert(etamin < etamax, "etamin must be smaller than etamax");
    float counter = 0;
    for (const auto& particle : mcParticles) {

      // primary
      if (selectPrimaries.value && !particle.isPhysicalPrimary()) {
        continue;
      }
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

  float countFT0A(const aod::McParticles& mcParticles) { return countMultInAcceptance<3.5f, 4.9f>(mcParticles); }
  float countFT0C(const aod::McParticles& mcParticles) { return countMultInAcceptance<-3.3f, -2.1f>(mcParticles); }
  float countFV0A(const aod::McParticles& mcParticles) { return countMultInAcceptance<2.2f, 5.1f>(mcParticles); }

  // Full tables (independent on central calibrations)
  void process(aod::McCollision const& mcCollision,
               aod::McParticles const& mcParticles)
  {
    const float nFT0M = countFT0A(mcParticles) + countFT0C(mcParticles);
    /*const float nFT0A = countFT0A(mcParticles);
    const float nFT0C = countFT0C(mcParticles);
    const float nFV0A = countFV0A(mcParticles);*/

    const float valueCentFT0M = h1dFT0M->GetBinContent(h1dFT0M->FindBin(nFT0M));
    /*const float valueCentFT0A = h1dFT0M->GetBinContent(h1dFT0M->FindBin(nFT0A));
    const float valueCentFT0C = h1dFT0M->GetBinContent(h1dFT0M->FindBin(nFT0C));
    const float valueCentFV0A = h1dFT0M->GetBinContent(h1dFT0M->FindBin(nFV0A));*/

    centFT0M(valueCentFT0M);
    /*centFT0A(valueCentFT0A);
    centFT0C(valueCentFT0C);
    centFV0A(valueCentFV0A);*/
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<mcCentrality>(cfgc)}; }
