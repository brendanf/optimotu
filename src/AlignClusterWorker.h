#ifndef OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED
#define OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED

#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>
#include "ClusterAlgorithm.h"

class AlignClusterWorker : public RcppParallel::Worker {
protected:
  const std::vector<std::string> &seq;
  ClusterAlgorithm &clust_algo;
  const uint8_t threads;
  std::mutex mutex;
  size_t _prealigned = 0, _aligned = 0;
  const bool verbose;
public :
  AlignClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    const uint8_t threads,
    bool verbose = false
  ) : seq(seq), clust_algo(clust_algo),
  threads(threads), verbose(verbose) {};

  size_t prealigned();

  size_t aligned();

  uint8_t nthreads();
};

#endif //OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED
