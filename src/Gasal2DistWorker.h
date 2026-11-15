// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_GASAL2DISTWORKER_H
#define OPTIMOTU_GASAL2DISTWORKER_H

#ifdef OPTIMOTU_R
#include "DistWorker.h"

class Gasal2DistWorker : public DistWorker {
protected:
  int match;
  int mismatch;
  int gap_open;
  int gap_extend;
public:
  Gasal2DistWorker(
    const std::vector<std::string> &seq,
    const double dist_threshold,
    DivisiblePairGenerator::Builder &pgb,
    SparseDistanceMatrix &sdm,
    int match,
    int mismatch,
    int gap_open,
    int gap_extend,
    int verbose = 0,
    std::size_t worker_threads = 0,
    Dissimilarity dissimilarity = Dissimilarity::IDENTITY
  );
  bool wants_gpu() const override;
  std::unique_ptr<GpuDistanceKernel> make_gpu_kernel() override;
};

template<int verbose>
class Gasal2DistWorkerImpl : public Gasal2DistWorker {
public:
  using Gasal2DistWorker::Gasal2DistWorker;
  void operator()(std::size_t begin, std::size_t end) override;
  void consume_gpu_tile(
    std::size_t pg_index,
    PairGenerator &pg,
    const GpuTileResult &result
  ) override;
};

std::unique_ptr<Gasal2DistWorker> create_gasal2_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder &pgb,
  SparseDistanceMatrix &sdm,
  int match,
  int mismatch,
  int gap_open,
  int gap_extend,
  int verbose = 0,
  std::size_t worker_threads = 0,
  Dissimilarity dissimilarity = Dissimilarity::IDENTITY
);

#endif
#endif
