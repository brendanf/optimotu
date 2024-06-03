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
    const uint8_t threads,
    const int match = 0, const int mismatch = 1,
    const int gap_open = 0, const int gap_extend = 1,
    const int gap_open2 = 0, const int gap_extend2 = 1,
    bool verbose = false
  );
};

class Wfa2SplitClusterWorker : public Wfa2ClusterWorker {
public:
  using Wfa2ClusterWorker::Wfa2ClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

class Wfa2ConcurrentClusterWorker : public Wfa2ClusterWorker {
public:
  using Wfa2ClusterWorker::Wfa2ClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

