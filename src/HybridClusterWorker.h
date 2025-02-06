#ifndef OPTIMOTU_HYBRIDCLUSTERWORKER_H
#define OPTIMOTU_HYBRIDCLUSTERWORKER_H

#include "AlignClusterWorker.h"

class HybridClusterWorker : public AlignClusterWorker {
protected:
  double breakpoint = 0.1;
public :
  HybridClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    const std::uint8_t threads,
    const double breakpoint = 0.1
  );
};

template<int verbose>
class HybridSplitClusterWorker : public HybridClusterWorker {
  using AlignClusterWorker::seq;
  using AlignClusterWorker::clust_algo;
  using AlignClusterWorker::threads;
  using AlignClusterWorker::mutex;
  using AlignClusterWorker::_prealigned;
  using AlignClusterWorker::_aligned;
public :
  using HybridClusterWorker::HybridClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

template<int verbose>
class HybridConcurrentClusterWorker : public HybridClusterWorker {
  using AlignClusterWorker::seq;
  using AlignClusterWorker::clust_algo;
  using AlignClusterWorker::threads;
  using AlignClusterWorker::mutex;
  using AlignClusterWorker::_prealigned;
  using AlignClusterWorker::_aligned;
public :
  using HybridClusterWorker::HybridClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

#endif
