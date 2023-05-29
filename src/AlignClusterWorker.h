#ifndef OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED
#define OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED

#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>
#include <bindings/cpp/WFAligner.hpp>
#include <edlib.h>
#include "pairwise_alignment.h"
#include "ClusterAlgorithm.h"
#include "ClusterMatrix.h"
#include "ClusterIndexedMatrix.h"
#include "ClusterTree.h"

class AlignClusterWorker : public RcppParallel::Worker {
protected:
  const std::vector<std::string> &seq;
  const double breakpoint;
  ClusterAlgorithm &clust_algo;
  const uint8_t threads;
  std::mutex mutex;
  size_t _prealigned = 0, _aligned = 0;
  const bool verbose;
public :
  AlignClusterWorker(
    const std::vector<std::string> &seq,
    const double breakpoint,
    ClusterAlgorithm &clust_algo,
    const uint8_t threads,
    bool verbose = false
  ) : seq(seq), breakpoint(breakpoint), clust_algo(clust_algo),
  threads(threads), verbose(verbose) {};

  size_t prealigned() {
    return _prealigned;
  }

  size_t aligned() {
    return _aligned;
  }
};

class HybridSplitClusterWorker : public AlignClusterWorker {
public :
  using AlignClusterWorker::AlignClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

class HybridConcurrentClusterWorker : public AlignClusterWorker {
public :
  using AlignClusterWorker::AlignClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

#endif //OPTIMOTU_ALIGNCLUSTERWORKER_H_INCLUDED
