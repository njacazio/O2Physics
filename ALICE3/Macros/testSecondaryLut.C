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

#include "ALICE3/Core/SecondaryLut.h"
#include "ALICE3/Core/SecondaryTrackSmearer.h"

#include <TH1.h>
#include <TSystem.h>
#include <TTree.h>
#include <TUUID.h>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

namespace
{
void requireSecondaryLut(bool condition, const char* message)
{
  if (!condition) {
    throw std::runtime_error(message);
  }
}
} // namespace

void testSecondaryLut()
{
  using Lut = o2::fastsim::SecondaryLut;
  Lut lut("SecondaryLut", 211, 0.139f, 20.f);
  requireSecondaryLut(lut.lookup(0, 0, 0, 0) == nullptr, "Uninitialized lookup");
  lut.setAxis(Lut::kNch, 2, 0.f, 2.f);
  lut.setAxis(Lut::kEta, 2, 0.f, 2.f);
  lut.setAxis(Lut::kPt, 3, 0.f, 3.f);
  lut.setAxis(Lut::kNHits, 4, -0.5f, 3.5f);
  lut.init();
  requireSecondaryLut(lut.size() == 48, "Unexpected number of bins");

  for (int nch = 0; nch < 2; ++nch) {
    for (int eta = 0; eta < 2; ++eta) {
      for (int hits = 0; hits < 4; ++hits) {
        for (int pt = 0; pt < 3; ++pt) {
          const auto idx = lut.index(nch, eta, pt, hits);
          requireSecondaryLut(idx == static_cast<std::size_t>(((nch * 2 + eta) * 4 + hits) * 3 + pt), "Incorrect index");
          auto& entry = lut.at(nch, eta, pt, hits);
          entry.valid = idx != 0;
          entry.eff = static_cast<float>(idx) / 100.f;
          entry.covm[0] = static_cast<float>(hits);
          requireSecondaryLut(lut.lookup(nch + 0.5f, eta + 0.5f, pt + 0.5f, hits) == &entry, "Incorrect lookup");
        }
      }
    }
  }
  requireSecondaryLut(lut.lookup(0.5f, 0.5f, 0.5f, -10.f) == &lut.at(0, 0, 0, 0), "Lower hit clamp");
  requireSecondaryLut(lut.lookup(0.5f, 0.5f, 0.5f, 99.f) == &lut.at(0, 0, 0, 3), "Upper hit clamp");

  o2::fastsim::SecondaryTrackSmearer smearer;
  requireSecondaryLut(smearer.adoptTable(211, lut), "Failed to adopt LUT");
  smearer.interpolateEfficiency(true);
  float interpolatedEff = -1.f;
  const auto* entry = smearer.getLUTEntry(211, 1.f, 1.5f, 2.5f, 2, interpolatedEff);
  requireSecondaryLut(entry == &smearer.getLUT(211)->at(1, 1, 2, 2) &&
                        std::abs(interpolatedEff - 0.32f) < 1.e-6f,
                      "Incorrect hit-dependent efficiency interpolation");
  requireSecondaryLut(std::abs(smearer.getEfficiency(211, 1.f, 1.5f, 2.5f, 3) - 0.35) < 1.e-6,
                      "Incorrect efficiency for different hit count");

  std::unique_ptr<TTree> tree(lut.toTree());
  requireSecondaryLut(tree->GetEntries() == 48 && tree->GetBranch("nHits") && tree->GetBranch("nHitsBin") &&
                        !tree->GetBranch("rad") && !tree->GetBranch("radBin") &&
                        !tree->GetBranch("prodR") && !tree->GetBranch("prodZ"),
                      "Incorrect tree structure");
  int nchBin = -1, etaBin = -1, ptBin = -1, hitsBin = -1;
  float hitsCenter = -1.f, efficiency = -1.f;
  tree->SetBranchAddress("nchBin", &nchBin);
  tree->SetBranchAddress("etaBin", &etaBin);
  tree->SetBranchAddress("ptBin", &ptBin);
  tree->SetBranchAddress("nHitsBin", &hitsBin);
  tree->SetBranchAddress("nHits", &hitsCenter);
  tree->SetBranchAddress("eff", &efficiency);
  for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
    tree->GetEntry(i);
    requireSecondaryLut(lut.index(nchBin, etaBin, ptBin, hitsBin) == static_cast<std::size_t>(i), "Incorrect tree bins");
    requireSecondaryLut(hitsCenter == hitsBin && efficiency == lut.entries()[i].eff, "Incorrect tree values");
  }

  std::unique_ptr<TH1D> projection(lut.projectEfficiency(Lut::kNHits));
  for (int hits = 0; hits < 4; ++hits) {
    double sum = 0.;
    int count = 0;
    for (int nch = 0; nch < 2; ++nch) {
      for (int eta = 0; eta < 2; ++eta) {
        for (int pt = 0; pt < 3; ++pt) {
          const auto& entry = lut.at(nch, eta, pt, hits);
          if (entry.valid) {
            sum += entry.eff;
            ++count;
          }
        }
      }
    }
    requireSecondaryLut(std::abs(projection->GetBinContent(hits + 1) - sum / count) < 1.e-6, "Incorrect hit projection");
  }

  const std::string fileName = std::string(gSystem->TempDirectory()) + "/secondary-lut-" + TUUID().AsString() + ".root";
  requireSecondaryLut(lut.writeToFile(fileName.c_str()), "Failed to write LUT");
  std::unique_ptr<Lut> restored(Lut::readFromFile(fileName.c_str()));
  gSystem->Unlink(fileName.c_str());
  requireSecondaryLut(restored && restored->size() == lut.size(), "Failed to read LUT");
  requireSecondaryLut(restored->getAxis(Lut::kNHits).nbins == 4 &&
                        restored->lookup(1.5f, 1.5f, 2.5f, 3.f)->covm[0] == 3.f,
                      "Hit axis or payload did not survive ROOT I/O");
}
