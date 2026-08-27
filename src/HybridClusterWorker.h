#ifndef OPTIMOTU_HYBRIDCLUSTERWORKER_H
#define OPTIMOTU_HYBRIDCLUSTERWORKER_H

#include "DistClusterWorker.h"

class HybridClusterWorker : public DistClusterWorker {
protected:
  double breakpoint = 0.1;
public :
  HybridClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    const double breakpoint = 0.1,
    int verbose = 0,
    std::size_t worker_threads = 0
  );
};

template<int verbose>
class HybridSplitClusterWorker : public HybridClusterWorker {
  using DistClusterWorker::seq;
  using DistClusterWorker::clust_algo;
  using DistClusterWorker::threads;
  using DistClusterWorker::mutex;
  using DistClusterWorker::_prealigned;
  using DistClusterWorker::_aligned;
  using DistClusterWorker::pair_generators;
public :
  HybridSplitClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    const double breakpoint = 0.1,
    std::size_t worker_threads = 0
  ) : HybridClusterWorker(seq, clust_algo, pgb, breakpoint, verbose, worker_threads) {};
  void operator()(std::size_t begin, std::size_t end);
};

template<int verbose>
class HybridConcurrentClusterWorker : public HybridClusterWorker {
  using DistClusterWorker::seq;
  using DistClusterWorker::clust_algo;
  using DistClusterWorker::threads;
  using DistClusterWorker::mutex;
  using DistClusterWorker::_prealigned;
  using DistClusterWorker::_aligned;
  using DistClusterWorker::pair_generators;
public :
  HybridConcurrentClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    const double breakpoint = 0.1,
    std::size_t worker_threads = 0
  ) : HybridClusterWorker(seq, clust_algo, pgb, breakpoint, verbose, worker_threads) {};
  void operator()(std::size_t begin, std::size_t end);
};

#endif
