// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED
#define OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED

#include <cstdint>
#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>
#include <algorithm>
#include "ClusterAlgorithm.h"
#include "PairGenerator.h"

class DistClusterWorker : public RcppParallel::Worker {
protected:
  const std::vector<std::string> &seq;
  ClusterAlgorithm &clust_algo;
  std::mutex mutex;
  size_t _prealigned = 0, _aligned = 0;
  std::vector<std::unique_ptr<PairGenerator>> pair_generators;
  const std::size_t threads;
public :
  DistClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pair_generator_builder,
    int verbose = 0,
    std::size_t worker_threads = 0
  ) : seq(seq), clust_algo(clust_algo),
  pair_generators(pair_generator_builder.build(verbose)),
  threads(
    worker_threads == 0
      ? std::max<std::size_t>(1, pair_generators.size())
      : std::max<std::size_t>(1, std::min(worker_threads, pair_generators.size()))
  ) {};

  size_t prealigned();

  size_t aligned();

  std::size_t n_threads();
};

#endif //OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED
