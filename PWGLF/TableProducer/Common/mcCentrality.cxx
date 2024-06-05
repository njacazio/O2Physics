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

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace o2::track;

/// Task to produce the response table
struct mcCentrality {
  using Trks = soa::Join<aod::TracksIU, aod::TracksExtra>;
  using Colls = aod::Collisions;

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

  HistogramRegistry histos{"histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(o2::framework::InitContext& initContext)
  {

    // Set up the CCDB
    const auto ts = ccdbTimestamp.value;
    ccdb->setURL(url.value);
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();
    ccdb->setCreatedNotAfter(std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count());
    if (ts >= 0) {
      LOGP(info, "Initialising LF TPC PID response for fixed timestamp {}:", ts);
      ccdb->setTimestamp(ts);
    }
    if (!fatalOnNonExisting) {
      LOG(warning) << "Setting the CCDB to not fatal when the object is not found. This means that you have no garantuee that the parameters used are correct. Expert user only!!!";
      ccdb->setFatalWhenNull(false);
    }
  }

  // Full tables (independent on central calibrations)
  void process(Colls const& collisions,
               Trks const& tracks,
               aod::BCsWithTimestamps const&)
  {
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc) { return WorkflowSpec{adaptAnalysisTask<mcCentrality>(cfgc)}; }
