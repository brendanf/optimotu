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
  using EdlibClusterWorker::EdlibClusterWorker;
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
  using EdlibClusterWorker::EdlibClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

#endif
