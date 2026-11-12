// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_AADISTWORKER_H_INCLUDED
#define OPTIMOTU_AADISTWORKER_H_INCLUDED

#ifdef OPTIMOTU_R

#include "DistWorker.h"
#include "PackedAaSequenceSet.h"
#include <memory>
#include <vector>

class AaDistWorker : public DistWorker {
protected:
  const PackedAaSequenceSet pss;
  const int min_overlap;
  const bool ignore_gaps;
  std::vector<float> lut_storage;
  const float *weight_lut;
  const bool return_all;

public:
  AaDistWorker(
    const std::vector<std::string> &seq,
    const double dist_threshold,
    DivisiblePairGenerator::Builder &pgb,
    SparseDistanceMatrix &sdm,
    const int min_overlap,
    const bool ignore_gaps,
    const float *weight_lut,
    bool return_all,
    int verbose = 0,
    std::size_t worker_threads = 0
  );
};

template<int verbose>
class AaDistWorkerImpl : public AaDistWorker {
  using DistWorker::dist_threshold;
  using DistWorker::pair_generators;
  using DistWorker::sdm;
  using AaDistWorker::pss;
  using AaDistWorker::min_overlap;
  using AaDistWorker::ignore_gaps;
  using AaDistWorker::weight_lut;
  using AaDistWorker::return_all;

public:
  using AaDistWorker::AaDistWorker;
  void operator()(std::size_t begin, std::size_t end) override;
};

std::unique_ptr<AaDistWorker> create_aa_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder &pgb,
  SparseDistanceMatrix &sdm,
  const int min_overlap,
  const bool ignore_gaps,
  const float *weight_lut,
  bool return_all,
  int verbose = 0,
  std::size_t worker_threads = 0
);

#endif
#endif
