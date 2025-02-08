// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED
#define OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED

#include <cstdint>
#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>
#include "ClusterAlgorithm.h"

class AlignClusterWorker : public RcppParallel::Worker {
protected:
  const std::vector<std::string> &seq;
  ClusterAlgorithm &clust_algo;
  const std::uint8_t threads;
  std::mutex mutex;
  size_t _prealigned = 0, _aligned = 0;
  const bool verbose;
public :
  AlignClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    const std::uint8_t threads,
    bool verbose = false
  ) : seq(seq), clust_algo(clust_algo),
  threads(threads), verbose(verbose) {};

  size_t prealigned();

  size_t aligned();

  std::uint8_t n_threads();
};

#endif //OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED
