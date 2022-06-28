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
/// \file   DetectorResponse.cxx
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>
/// \since  2020-07-30
/// \brief  Handler for any detector (or other entity) response.
///         This provides the basic quantities computed by any response i.e. expected values, resolutions and Nsigmas
///

// ROOT includes
#include "Rtypes.h"
#include "TMath.h"
#include "TFile.h"

// O2 includes
#include "Framework/Logger.h"
#include "ReconstructionDataFormats/PID.h"
#include "PID/ParamBase.h"
#include "DetectorResponse.h"

using namespace o2::pid;

void DetectorResponse::LoadParamFromCCDB(const Param_t ptype, const std::string& path, long& timestamp)
{
  mParamCCDBPaths[ptype] = path;
  LoadParam(ptype, mCCDBManager->getForTimeStamp<Parametrization>(path, timestamp));
}

void DetectorResponse::UpdateParamFromCCDB(const Param_t ptype, long& timestamp)
{
  for (Param_t i = kSignal; i < kNParams; i++) {
    if (!mParam[i]) {
      continue;
    }
    LoadParamFromCCDB(ptype, mParamCCDBPaths[ptype], timestamp);
  }
}

void DetectorResponse::LoadParam(const Param_t ptype, Parametrization* param)
{
  if (!param) {
    LOG(fatal) << "Got no param for parametrization " << ParamName[ptype];
  }
  mParam[ptype] = param;
}

void DetectorResponse::LoadParamFromFile(const std::string& fname, const std::string& pname, const Param_t ptype)
{
  TFile f(fname.c_str(), "READ");
  if (!f.Get(pname.c_str())) {
    LOG(fatal) << "Did not find parametrization " << pname << " in file " << fname;
  }
  LOG(info) << "Loading parametrization " << pname << " from TFile " << fname;
  f.GetObject(pname.c_str(), mParam[ptype]);
  f.Close();
  mParam[ptype]->Print();
}

void DetectorResponse::SetParameters(const DetectorResponse::Param_t ptype, std::vector<pidvar_t> p)
{
  if (!mParam[ptype]) {
    const std::string pname = std::string(ParamName[ptype]) + "_default_param";
    LOG(info) << "Creating new parametrization " << pname << " of size " << p.size();
    mParam[ptype] = new Parametrization(pname, p.size());
    mParam[ptype]->Print();
  }
  mParam[ptype]->SetParameters(p);
}
