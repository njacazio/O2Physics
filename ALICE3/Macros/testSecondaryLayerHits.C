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

#include "ALICE3/Core/FastTracker.h"

#include <CommonConstants/MathConstants.h>

#include <array>
#include <cmath>
#include <stdexcept>

void testSecondaryLayerHits()
{
  using o2::fastsim::DetLayer;
  o2::fastsim::FastTracker tracker;
  tracker.addLayer("passive", 1.f, 100.f, 0.f, 0.f, 0.f, 0.f, 0.f, DetLayer::kLayerInert);
  tracker.addLayer("inner", 2.f, 3.f, 0.f, 0.f, 0.f, 0.f, 1.f, DetLayer::kLayerSilicon);
  tracker.addLayer("middle", 6.f, 4.f, 0.f, 0.f, 0.f, 0.f, 1.f, DetLayer::kLayerSilicon);
  tracker.addLayer("outer", 11.f, 20.f, 0.f, 0.f, 0.f, 0.f, 1.f, DetLayer::kLayerSilicon);
  tracker.addLayer("tof", 12.f, 100.f, 0.f, 0.f, 0.f, 0.f, 1.f, DetLayer::kLayerTOF);
  tracker.setMagneticField(0.f);

  const o2::track::TrackParCov::covMat_t covariance{};
  const o2::track::TrackParCov::params_t inwardParams{0.f, 0.f, 0.f, 0.5f, 1.f};
  o2::track::TrackParCov inward(-4.f, o2::constants::math::PI, inwardParams, covariance);
  auto hits = tracker.countSiliconLayerHits(inward);
  if (hits.count != 2 || hits.innerRadius != 2.f || hits.outerRadius != 11.f) {
    throw std::runtime_error("Incorrect displaced inward straight-track intersections");
  }

  tracker.setMagneticField(20.f);
  const float qOverPt = 1.f / (5.f * 20.f * o2::constants::math::B2C);
  for (const float sign : {-1.f, 1.f}) {
    const o2::track::TrackParCov::params_t curlingParams{0.f, 0.f, 0.f, 0.f, sign * qOverPt};
    o2::track::TrackParCov curling(0.f, 0.f, curlingParams, covariance);
    hits = tracker.countSiliconLayerHits(curling);
    if (hits.count != 2 || hits.innerRadius != 2.f || hits.outerRadius != 6.f) {
      throw std::runtime_error("Incorrect helical silicon-layer intersections");
    }
    const o2::track::TrackParCov::params_t forwardParams{0.f, 0.f, 0.f, 1.f, sign * qOverPt};
    o2::track::TrackParCov forward(0.f, 0.f, forwardParams, covariance);
    hits = tracker.countSiliconLayerHits(forward);
    if (hits.count != 1 || hits.innerRadius != 2.f || hits.outerRadius != 2.f) {
      throw std::runtime_error("Incorrect helical longitudinal acceptance");
    }
  }
}
