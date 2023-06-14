#include "ClusterWorker.h"

ClusterWorker::ClusterWorker(std::istream &file, const int threads) :
  file(file), threads(threads) {}

int ClusterWorker::n_threads() {
  return this->threads;
}

MergeClusterWorker::MergeClusterWorker(
  ClusterAlgorithm *algo,
  std::istream &file,
  const int threads
) :
  ClusterWorker(file, threads),
  algo_list(threads) {
  // OPTIMOTU_COUT << "MergeClusterWorker constructor start...";
  for (int i = 0; i < threads; ++i) {
    algo_list[i] = algo->make_child();
  }
  // OPTIMOTU_COUT << "done" << std::endl;
}

void MergeClusterWorker::operator()(size_t begin, size_t end) {
  DistanceElement d;
  std::vector<DistanceElement> buffer;
  buffer.reserve(100);
  // OPTIMOTU_COUT << "Starting MergeClusterWorker thread " << begin << std::endl;
  while (true) {
    mutex.lock();
    for (int i = 0; i < 100 && file; ++i) {
      file >> d;
      buffer.push_back(d);
    }
    mutex.unlock();
    if (buffer.size() == 0) break;
    for (auto d : buffer) {
      algo_list[begin]->operator()(d, begin);
    }
    buffer.clear();
  }
  // mutex.lock();
  // OPTIMOTU_COUT << "ClusterWorker thread " << begin << " merging..." << std::endl;
  // mutex.unlock();
  algo_list[begin]->merge_into_parent();

  // mutex.lock();
  // OPTIMOTU_COUT << "ClusterWorker thread " << begin << " exiting" << std::endl;
  // mutex.unlock();
}

void MergeClusterWorker::finalize() {}

ConcurrentClusterWorker::ConcurrentClusterWorker(
  ClusterAlgorithm *algo,
  std::istream &file,
  const int threads
) :
  ClusterWorker(file, threads),
  algo(algo) {}

void ConcurrentClusterWorker::operator()(size_t begin, size_t end) {
  DistanceElement d;
  std::vector<DistanceElement> buffer;
  buffer.reserve(100);
  // OPTIMOTU_COUT << "Starting ConcurrentClusterWorker thread " << begin << std::endl;
  while (true) {
    mutex.lock();
    // OPTIMOTU_COUT << "buffering..." << std::flush;
    for (int i = 0; i < 100 && file; ++i) {
      file >> d;
      buffer.push_back(d);
    }
    // OPTIMOTU_COUT << "done" << std::endl;
    mutex.unlock();
    if (buffer.size() == 0) break;
    for (auto d : buffer) {
      algo->operator()(d, begin);
    }
    buffer.clear();
  }

  // mutex.lock();
  // OPTIMOTU_COUT << "ConcurrentClusterWorker thread " << begin << " exiting" << std::endl;
  // mutex.unlock();
}

void ConcurrentClusterWorker::finalize() {}

HierarchicalClusterWorker::HierarchicalClusterWorker(
  ClusterAlgorithm *algo,
  std::istream &file,
  const int threads,
  const int shards
) : ClusterWorker(file, threads), algo_list(shards),
thread_count(shards), shards(shards) {
  // OPTIMOTU_COUT << "HieararchicalClusterWorker constructor start...";
  for (int i = 0; i < shards; ++i) {
    algo_list[i] = algo->make_child();
  }
  // OPTIMOTU_COUT << "done" << std::endl;
}

void HierarchicalClusterWorker::operator()(size_t begin, size_t end) {
  DistanceElement d;
  std::vector<DistanceElement> buffer;
  buffer.reserve(100);
  size_t i = begin % shards;
  // mutex.lock();
  // size_t tc =
    ++thread_count[i];
  // OPTIMOTU_COUT << "Starting HierarchicalClusterWorker thread " << begin
  //           << " (shard " << i
  //           << " with " << tc
  //           << " previous threads)"
  //           << std::endl;
  // mutex.unlock();
  while (true) {
    mutex.lock();
    for (int i = 0; i < 100 && file; ++i) {
      file >> d;
      buffer.push_back(d);
    }
    mutex.unlock();
    if (buffer.size() == 0) break;
    for (auto d : buffer) {
      algo_list[i]->operator()(d, begin);
    }
    buffer.clear();
  }
  if (--thread_count[i] == 0) {
    // mutex.lock();
    // OPTIMOTU_COUT << "HierarchicalClusterWorker thread " << begin
    //           << " merging shard " << i
    //           << "..." << std::endl;
    // mutex.unlock();
    algo_list[i]->merge_into_parent();
  }

    // mutex.lock();
    // OPTIMOTU_COUT << "HieararchicalClusterWorker thread " << begin
    //           << " exiting" << std::endl;
    // mutex.unlock();
  }

  void HierarchicalClusterWorker::finalize() {

}
