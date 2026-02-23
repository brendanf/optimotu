// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED
#define OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED

#include <cstdint>
#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>
#include "ClusterAlgorithm.h"
#include "PairGenerator.h"

class DistClusterWorker : public RcppParallel::Worker {
protected:
  const std::vector<std::string> &seq;
  ClusterAlgorithm &clust_algo;
  std::mutex mutex;
  size_t _prealigned = 0, _aligned = 0;
  std::vector<std::unique_ptr<PairGenerator>> pair_generators;
  const std::uint8_t threads;
public :
  DistClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pair_generator_builder
  ) : seq(seq), clust_algo(clust_algo),
  pair_generators(pair_generator_builder.build()),
  threads(pair_generators.size()) {};

  size_t prealigned();

  size_t aligned();

  std::uint8_t n_threads();
};

#endif //OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED
