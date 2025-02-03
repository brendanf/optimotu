#include "AlignClusterWorker.h"

class Wfa2ClusterWorker : public AlignClusterWorker {
protected:
  int match = 0, mismatch = 1,
    gap_open = 0, gap_extend = 1,
    gap_open2 = 0, gap_extend2 = 1;
public :
  Wfa2ClusterWorker(
    const std::vector<std::string> &seq,
    ClusterAlgorithm &clust_algo,
    const std::uint8_t threads,
    const int match = 0, const int mismatch = 1,
    const int gap_open = 0, const int gap_extend = 1,
    const int gap_open2 = 0, const int gap_extend2 = 1
  );
};

template<int verbose>
class Wfa2SplitClusterWorker : public Wfa2ClusterWorker {
  using AlignClusterWorker::seq;
  using AlignClusterWorker::clust_algo;
  using AlignClusterWorker::threads;
  using AlignClusterWorker::mutex;
  using AlignClusterWorker::_prealigned;
  using AlignClusterWorker::_aligned;

  using Wfa2ClusterWorker::match;
  using Wfa2ClusterWorker::mismatch;
  using Wfa2ClusterWorker::gap_open;
  using Wfa2ClusterWorker::gap_extend;
  using Wfa2ClusterWorker::gap_open2;
  using Wfa2ClusterWorker::gap_extend2;
public:
  using Wfa2ClusterWorker::Wfa2ClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

template<int verbose>
class Wfa2ConcurrentClusterWorker : public Wfa2ClusterWorker {
  using AlignClusterWorker::seq;
  using AlignClusterWorker::clust_algo;
  using AlignClusterWorker::threads;
  using AlignClusterWorker::mutex;
  using AlignClusterWorker::_prealigned;
  using AlignClusterWorker::_aligned;

  using Wfa2ClusterWorker::match;
  using Wfa2ClusterWorker::mismatch;
  using Wfa2ClusterWorker::gap_open;
  using Wfa2ClusterWorker::gap_extend;
  using Wfa2ClusterWorker::gap_open2;
  using Wfa2ClusterWorker::gap_extend2;
public:
  using Wfa2ClusterWorker::Wfa2ClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

