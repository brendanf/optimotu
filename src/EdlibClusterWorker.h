#include "AlignClusterWorker.h"

class EdlibClusterWorker : public AlignClusterWorker {
public :
  using AlignClusterWorker::AlignClusterWorker;
};

template <int verbose>
class EdlibSplitClusterWorker : public EdlibClusterWorker {
  using AlignClusterWorker::seq;
  using AlignClusterWorker::clust_algo;
  using AlignClusterWorker::threads;
  using AlignClusterWorker::mutex;
  using AlignClusterWorker::_prealigned;
  using AlignClusterWorker::_aligned;
public:
  using EdlibClusterWorker::EdlibClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

template <int verbose>
class EdlibConcurrentClusterWorker : public EdlibClusterWorker {
  using AlignClusterWorker::seq;
  using AlignClusterWorker::threads;
  using AlignClusterWorker::mutex;
  using AlignClusterWorker::clust_algo;
  using AlignClusterWorker::_prealigned;
  using AlignClusterWorker::_aligned;
public:
  using EdlibClusterWorker::EdlibClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};
