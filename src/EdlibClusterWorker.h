#include "AlignClusterWorker.h"

class EdlibClusterWorker : public AlignClusterWorker {
public :
  using AlignClusterWorker::AlignClusterWorker;
};

class EdlibSplitClusterWorker : public EdlibClusterWorker {
public:
  using EdlibClusterWorker::EdlibClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};

class EdlibConcurrentClusterWorker : public EdlibClusterWorker {
public:
  using EdlibClusterWorker::EdlibClusterWorker;
  void operator()(std::size_t begin, std::size_t end);
};
