#ifndef OPTIMOTU_EDLIBCLUSTERWORKER_H
#define OPTIMOTU_EDLIBCLUSTERWORKER_H

#include "DistClusterWorker.h"

class EdlibClusterWorker : public DistClusterWorker {
public :
  using DistClusterWorker::DistClusterWorker;
};

template <int verbose>
class EdlibSplitClusterWorker : public EdlibClusterWorker {
  using DistClusterWorker::seq;
  using DistClusterWorker::clust_algo;
  using DistClusterWorker::threads;
  using DistClusterWorker::mutex;
  using DistClusterWorker::_prealigned;
  using DistClusterWorker::_aligned;
public:
  EdlibSplitClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    std::size_t worker_threads = 0
  ) : EdlibClusterWorker(seq, clust_algo, pgb, verbose, worker_threads) {};
  void operator()(std::size_t begin, std::size_t end);
};

template <int verbose>
class EdlibConcurrentClusterWorker : public EdlibClusterWorker {
  using DistClusterWorker::seq;
  using DistClusterWorker::threads;
  using DistClusterWorker::mutex;
  using DistClusterWorker::clust_algo;
  using DistClusterWorker::_prealigned;
  using DistClusterWorker::_aligned;
  public:
  EdlibConcurrentClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    DivisiblePairGenerator::Builder & pgb,
    std::size_t worker_threads = 0
  ) : EdlibClusterWorker(seq, clust_algo, pgb, verbose, worker_threads) {};
  void operator()(std::size_t begin, std::size_t end);
};

#endif
