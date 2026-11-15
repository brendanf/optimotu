// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_GASAL2CLUSTERWORKER_H
#define OPTIMOTU_GASAL2CLUSTERWORKER_H

#include "DistClusterWorker.h"
#include "gpu_tile.h"

class Gasal2ClusterWorker : public DistClusterWorker {
protected:
  int match = 0;
  int mismatch = 1;
  int gap_open = 0;
  int gap_extend = 1;
public:
  Gasal2ClusterWorker(
    const SequenceSet &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder &pgb,
    int match,
    int mismatch,
    int gap_open,
    int gap_extend,
    int verbose,
    std::size_t worker_threads,
    Dissimilarity dissimilarity
  );
  bool wants_gpu() const override;
  std::unique_ptr<GpuDistanceKernel> make_gpu_kernel() override;
};

template<int verbose>
class Gasal2SplitClusterWorker : public Gasal2ClusterWorker {
public:
  using Gasal2ClusterWorker::Gasal2ClusterWorker;
  void operator()(std::size_t begin, std::size_t end) override;
  void consume_gpu_tile(
    std::size_t pg_index,
    PairGenerator &pg,
    const GpuTileResult &result
  ) override;
};

template<int verbose>
class Gasal2ConcurrentClusterWorker : public Gasal2ClusterWorker {
public:
  using Gasal2ClusterWorker::Gasal2ClusterWorker;
  void operator()(std::size_t begin, std::size_t end) override;
  void consume_gpu_tile(
    std::size_t pg_index,
    PairGenerator &pg,
    const GpuTileResult &result
  ) override;
};

#endif
