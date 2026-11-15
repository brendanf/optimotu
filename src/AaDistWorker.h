// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_AADISTWORKER_H_INCLUDED
#define OPTIMOTU_AADISTWORKER_H_INCLUDED

#ifdef OPTIMOTU_R

#include "DistWorker.h"
#include <tuple>
#include "EvoModel.h"
#include "PackedAaSequenceSet.h"
#include "gpu_tile.h"
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
  bool has_evo = false;
  EvoModel evo;

  void score_targets(
    int query,
    const uint8_t *pos_major,
    const int *target_idx,
    int n_target,
    int n_cols,
    double dist_threshold,
    bool early_exit,
    double *dist_out,
    int *n_ok_out
  ) const;
  std::tuple<bool, int, double> score_pair(int i, int j) const;
  bool wants_gpu() const override;
  std::unique_ptr<GpuDistanceKernel> make_gpu_kernel() override;

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
    std::size_t worker_threads = 0,
    const EvoModel *evo_model = nullptr
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
  using AaDistWorker::score_targets;
  using AaDistWorker::score_pair;

public:
  using AaDistWorker::AaDistWorker;
  void operator()(std::size_t begin, std::size_t end) override;
  void consume_gpu_tile(
    std::size_t pg_index,
    PairGenerator &pg,
    const GpuTileResult &result
  ) override;
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
  std::size_t worker_threads = 0,
  const EvoModel *evo_model = nullptr
);

#endif
#endif
