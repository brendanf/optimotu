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
public :
  AlignClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    const std::uint8_t threads
  ) : seq(seq), clust_algo(clust_algo),
  threads(threads) {};

  size_t prealigned();

  size_t aligned();

  std::uint8_t n_threads();
};

#endif //OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED
