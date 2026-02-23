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
    const double breakpoint = 0.1
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
  using HybridClusterWorker::HybridClusterWorker;
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
  using HybridClusterWorker::HybridClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

#endif
