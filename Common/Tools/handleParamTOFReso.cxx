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
/// \file   handleParamTOFReso.cxx
/// \author Nicolò Jacazio nicolo.jacazio@cern.ch
/// \since  2020-06-22
/// \brief  A simple tool to produce Bethe Bloch parametrization objects for the TOF PID Response
///

#include "TGraph.h"
#include "TSystem.h"
#include "TCanvas.h"

#include "Common/Core/PID/TOFReso.h"
#include "Common/Core/PID/DetectorResponse.h"
#include "Common/Core/PID/PIDTOF.h"
#include "handleParamBase.h"

#include <chrono>
using namespace std::chrono;
using namespace o2::pid::tof;
using namespace o2::pid;

struct DebugTrack { // Track that mimics the O2 data structure
  float mp = 0.1f;
  float p() const { return mp; }
  bool isEvTimeDefined() const { return true; }
  bool hasTOF() const { return true; }
  float tofEvTime() const { return 0.f; }
  float tofEvTimeErr() const { return 20.f; }
  float length() const { return 400.f; }
  float tofSignal() const { return length() * 1.01 / kCSPEED; }
  int trackType() const { return 0; };
  float tofExpMom() const
  {
    constexpr float mass = o2::constants::physics::MassPionCharged; // default pid = pion
    const float expBeta = (length() / (tofSignal() * kCSPEED));
    return mass * expBeta / std::sqrt(1.f - expBeta * expBeta);
  }

} debugTrack;

bool initOptionsAndParse(bpo::options_description& options, int argc, char* argv[])
{
  options.add_options()(
    "url,u", bpo::value<std::string>()->default_value("http://alice-ccdb.cern.ch"), "URL of the CCDB database e.g. http://ccdb-test.cern.ch:8080 or http://alice-ccdb.cern.ch")(
    "ccdb-path,c", bpo::value<std::string>()->default_value("Analysis/PID/TOF"), "CCDB path for storage/retrieval")(
    "rct-path", bpo::value<std::string>()->default_value("RCT/Info/RunInformation"), "path to the ccdb RCT objects for the SOR/EOR timestamps")(
    "start,s", bpo::value<long>()->default_value(0), "Start timestamp of object validity. If 0 and runnumber != 0 it will be set to the run SOR")(
    "stop,S", bpo::value<long>()->default_value(0), "Stop timestamp of object validity. If 0 and runnumber != 0 it will be set to the run EOR")(
    "timestamp,T", bpo::value<long>()->default_value(-1), "Timestamp of the object to retrieve, used in alternative to the run number")(
    "runnumber,R", bpo::value<unsigned int>()->default_value(0), "Timestamp of the object to retrieve, used in alternative to the timestamp (if 0 using the timestamp)")(
    "delete-previous,delete_previous,d", bpo::value<int>()->default_value(0), "Flag to delete previous versions of converter objects in the CCDB before uploading the new one so as to avoid proliferation on CCDB")(
    "save-to-file,file,f,o", bpo::value<std::string>()->default_value(""), "Option to save parametrization to file instead of uploading to ccdb")(
    "read-from-file,i", bpo::value<std::string>()->default_value(""), "Option to get parametrization from a file")(
    "reso-name,n", bpo::value<std::string>()->default_value("TOFReso"), "Name of the parametrization object")(
    "mode,m", bpo::value<unsigned int>()->default_value(1), "Working mode: 0 push 1 pull and test")(
    "p0", bpo::value<float>()->default_value(0.008f), "Parameter 0 of the TOF resolution")(
    "p1", bpo::value<float>()->default_value(0.008f), "Parameter 1 of the TOF resolution")(
    "p2", bpo::value<float>()->default_value(0.002f), "Parameter 2 of the TOF resolution")(
    "p3", bpo::value<float>()->default_value(40.0f), "Parameter 3 of the TOF resolution")(
    "p4", bpo::value<float>()->default_value(60.0f), "Parameter 4 of the TOF resolution: average TOF resolution")(
    "verbose,v", bpo::value<int>()->default_value(0), "Verbose level 0, 1")(
    "help,h", "Produce help message.");
  try {
    bpo::store(parse_command_line(argc, argv, options), arguments);

    // help
    if (arguments.count("help")) {
      LOG(info) << options;
      return false;
    }

    bpo::notify(arguments);
  } catch (const bpo::error& e) {
    LOG(error) << e.what() << "\n";
    LOG(error) << "Error parsing command line arguments; Available options:";
    LOG(error) << options;
    return false;
  }
  return true;
}

int main(int argc, char* argv[])
{
  bpo::options_description options("Allowed options");
  if (!initOptionsAndParse(options, argc, argv)) {
    return 1;
  }

  // Fetch options
  const auto mode = arguments["mode"].as<unsigned int>();
  const auto runnumber = arguments["runnumber"].as<unsigned int>();
  auto timestamp = arguments["timestamp"].as<long>();
  const auto path = arguments["ccdb-path"].as<std::string>();
  auto start = arguments["start"].as<long>();
  auto stop = arguments["stop"].as<long>();

  // Init CCDB
  const std::string url = arguments["url"].as<std::string>();
  api.init(url);
  if (!api.isHostReachable()) {
    LOG(warning) << "CCDB host " << url << " is not reacheable, cannot go forward";
    return 1;
  }

  // Init timestamps
  setupTimestamps(timestamp, start, stop);

  TOFReso* resoOld = nullptr;
  TOFResoParams* reso = nullptr;
  const std::string reso_name = arguments["reso-name"].as<std::string>();
  if (mode == 0) { // Push mode
    LOG(info) << "Handling TOF parametrization in create mode";
    const std::string input_file_name = arguments["read-from-file"].as<std::string>();
    reso = new TOFResoParams();
    if (!input_file_name.empty()) { // Load parameters from input file
      LOG(info) << "Loading parameters from file";
      reso->LoadParamFromFile(input_file_name.c_str(), reso_name.c_str());
    } else { // Create new object
      LOG(info) << "Loading parameters from command line";
      reso->SetParameters(std::array<float, 5>{arguments["p0"].as<float>(),
                                               arguments["p1"].as<float>(),
                                               arguments["p2"].as<float>(),
                                               arguments["p3"].as<float>(),
                                               arguments["p4"].as<float>()});
    }
    reso->Print();
    const std::string fname = arguments["save-to-file"].as<std::string>();
    if (!fname.empty()) { // Saving it to file
      LOG(info) << "Saving parametrization to file " << fname;
      TFile f(fname.data(), "RECREATE");
      reso->Write();
      f.ls();
      f.Close();
    } else { // Saving it to CCDB
      LOG(info) << "Saving parametrization to CCDB " << path << " with validity " << start << " "
                << stop;

      std::map<std::string, std::string> metadata;
      if (runnumber != 0) {
        metadata["runnumber"] = Form("%i", runnumber);
      }
      reso->AddToMetadata(metadata);
      // Storing parametrization parameters
      storeOnCCDB(path + "/Parameters/" + reso_name, metadata, start, stop, reso);
    }
  } else if (mode == 1) { // Pull and test mode
    LOG(info) << "Handling TOF parametrization in test mode for timestamp "
              << timestamp << " -> " << timeStampToHReadble(timestamp);
    reso = retrieveFromCCDB<TOFResoParams>(path + "/" + reso_name, timestamp);
    reso->Print();
    using RespImp = ExpTimes<DebugTrack, 2>;
    LOG(info) << "TOF expected resolution at p=" << debugTrack.p() << " GeV/c and mass " << RespImp::mMassZ << ":" << RespImp::GetExpectedSigma(*reso, debugTrack);
  } else { // Create and test + performance
    LOG(info) << "Creating TOF parametrization and testing";
    reso = new TOFResoParams();
    reso->Print();
    DetectorResponse response;
    if (!resoOld) {
      resoOld = new TOFReso();
      const std::vector<float> resoparams = {arguments["p0"].as<float>(), arguments["p1"].as<float>(), arguments["p2"].as<float>(), arguments["p3"].as<float>(), arguments["p4"].as<float>()};
      resoOld->SetParameters(resoparams);
    }
    response.LoadParam(DetectorResponse::kSigma, resoOld);
    // Draw it
    using RespImp = ExpTimes<DebugTrack, 2>;
    //
    std::map<std::string, TGraph*> graphs;
    graphs["ExpSigma"] = new TGraph();
    graphs["NSigma"] = new TGraph();
    graphs["durationExpSigma"] = new TGraph();
    graphs["durationNSigma"] = new TGraph();
    //
    graphs["ExpSigmaOld"] = new TGraph();
    graphs["NSigmaOld"] = new TGraph();
    graphs["durationExpSigmaOld"] = new TGraph();
    graphs["durationNSigmaOld"] = new TGraph();
    const int nsamp = 1000;
    for (int i = 0; i < nsamp; i++) {
      debugTrack.mp += 0.01f;
      //
      auto start = high_resolution_clock::now();
      RespImp::GetExpectedSigma(*reso, debugTrack);
      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<nanoseconds>(stop - start).count();
      graphs["ExpSigma"]->SetPoint(i, debugTrack.p(), RespImp::GetExpectedSigma(*reso, debugTrack));
      graphs["durationExpSigma"]->SetPoint(i + 1, i, duration);
      //
      start = high_resolution_clock::now();
      RespImp::GetSeparation(*reso, debugTrack);
      stop = high_resolution_clock::now();
      duration = duration_cast<nanoseconds>(stop - start).count();
      graphs["NSigma"]->SetPoint(i, debugTrack.p(), RespImp::GetSeparation(*reso, debugTrack));
      graphs["durationNSigma"]->SetPoint(i + 1, i, duration);
      //
      start = high_resolution_clock::now();
      RespImp::GetExpectedSigma(response, debugTrack);
      stop = high_resolution_clock::now();
      duration = duration_cast<nanoseconds>(stop - start).count();
      graphs["ExpSigmaOld"]->SetPoint(i, debugTrack.p(), RespImp::GetExpectedSigma(*reso, debugTrack));
      graphs["durationExpSigmaOld"]->SetPoint(i + 1, i, duration);
      //
      start = high_resolution_clock::now();
      RespImp::GetSeparation(response, debugTrack);
      stop = high_resolution_clock::now();
      duration = duration_cast<nanoseconds>(stop - start).count();
      graphs["NSigmaOld"]->SetPoint(i, debugTrack.p(), RespImp::GetSeparation(*reso, debugTrack));
      graphs["durationNSigmaOld"]->SetPoint(i + 1, i, duration);
    }
    TFile fdebug("/tmp/tofParamDebug.root", "UPDATE");
    TString dn = Form("%i", fdebug.GetListOfKeys().GetEntries());
    fdebug.mkdir(dn);
    fdebug.cd(dn);
    for (const auto& i : graphs) {
      i.second->SetName(i.first.c_str());
      i.second->SetTitle(i.first.c_str());
      i.second->Write();
    }
    fdebug.Close();
    LOG(info) << "TOF expected resolution at p=" << debugTrack.p() << " GeV/c and mass " << RespImp::mMassZ << ": " << RespImp::GetExpectedSigma(*reso, debugTrack) << " ps";
  }

  return 0;
}
