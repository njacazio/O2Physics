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
/// \file   DetectorResponse.h
/// \author Nicolò Jacazio <nicolo.jacazio@cern.ch>
/// \since  2020-07-30
/// \brief  Handler for any detector (or other entity) response.
///         This provides the basic quantities computed by any response i.e. expected values, resolutions and Nsigmas
///

#ifndef O2_ANALYSIS_PID_DETECTORRESPONSE_H_
#define O2_ANALYSIS_PID_DETECTORRESPONSE_H_

#include <Rtypes.h>
#include <array>
#include <vector>
#include "Common/Core/PID/ParamBase.h"
#include <CCDB/BasicCCDBManager.h>

namespace o2::pid
{
/// \brief Class to handle the general detector response
class DetectorResponse
{
 public:
  DetectorResponse() = default;
  virtual ~DetectorResponse() = default;

  /// Enumeration of the different types of parametrizations

  // enum class Param_t : uint8_t {
  //   kSignal = 0,
  //   kSigma,
  //   kNParams
  // };
  typedef uint8_t Param_t;
  static constexpr Param_t kSignal = 0;
  static constexpr Param_t kSigma = 1;
  static constexpr Param_t kNParams = 2;

  static constexpr std::array<char const*, kNParams> ParamName = {{"Signal", "Sigma"}};

  /// Setter for the parametrization from input TFile
  /// \param fname File name used for input
  /// \param pname Name of the parametrization in the file
  /// \param ptype Type of the parametrization
  void LoadParamFromFile(const std::string& fname, const std::string& pname, const Param_t ptype);

  /// Setter for the parametrization
  /// \param ptype Type of the parametrization
  /// \param param Parametrization
  void LoadParam(const Param_t ptype, Parametrization* param);

  /// Setter for the parametrization from the CCDB
  /// \param ptype Type of the parametrization
  /// \param path Path to the parametrization object
  /// \param timestamp timestamp of the parametrization object
  void LoadParamFromCCDB(const Param_t ptype, const std::string& path, long& timestamp);

  /// Setter for the parametrization from the CCDB, in update mode, only if object is not valid anymore
  /// \param ptype Type of the parametrization
  /// \param timestamp timestamp of the parametrization object
  void UpdateParamFromCCDB(const Param_t ptype, long& timestamp);

  /// Setter for the CCDB handler to use for the loading of the parametrization
  /// \param ccdb Configured CCDB object
  void InitCCDBManager(o2::ccdb::BasicCCDBManager& ccdb) { mCCDBManager = &ccdb; }

  /// Getter for the parametrizations
  Parametrization* GetParam(const Param_t ptype) const { return mParam[ptype]; }

  /// Setter for the parametrizations parameters, if the parametrization is not yet initialized a new parametrization is created without any implementation and just parameters
  /// \param ptype parametrization type
  /// \param p vector with parameters
  void SetParameters(const Param_t ptype, std::vector<pidvar_t> p);

  /// Getter for the value of the parametrization
  /// \param ptype parametrization type
  /// \param x array with parameters
  pidvar_t operator()(const Param_t ptype, const pidvar_t* x) const { return mParam[ptype]->operator()(x); }

 private:
  std::array<Parametrization*, kNParams> mParam; /// Parametrizations for the expected signal and sigma
  // Data members for auto updating the parametrizations
  std::array<std::string, kNParams> mParamCCDBPaths; /// Array of paths on the CCDB used to fetch parametrizations, used in update mode
  o2::ccdb::BasicCCDBManager* mCCDBManager;          /// Pointer to a configured CCDB manager, used to fetch and update parametrizations
};

} // namespace o2::pid

#endif // O2_ANALYSIS_PID_DETECTORRESPONSE_H_
