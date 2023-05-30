#include "ClusterWorker.h"

ClusterWorker::ClusterWorker(std::istream &file, const int threads) : file(file), threads(threads) {};

int ClusterWorker::n_threads() {
  return this->threads;
}

MergeClusterWorker::MergeClusterWorker(ClusterAlgorithm *algo, std::istream &file, const int threads) :
  ClusterWorker(file, threads),
  algo_list(threads) {
  std::cout << "MergeClusterWorker constructor start...";
  for (int i = 0; i < threads; ++i) {
    algo_list[i] = algo->make_child();
  }
  std::cout << "done" << std::endl;
};

void MergeClusterWorker::operator()(size_t begin, size_t end) {
  DistanceElement d;
  std::vector<DistanceElement> buffer;
  buffer.reserve(100);
  std::cout << "Starting ClusterWorker thread " << begin << std::endl;
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
  mutex.lock();
  std::cout << "ClusterWorker thread " << begin << " merging..." << std::endl;
  mutex.unlock();
  algo_list[begin]->merge_into_parent();

  mutex.lock();
  std::cout << "ClusterWorker thread " << begin << " exiting" << std::endl;
  mutex.unlock();
}

void MergeClusterWorker::finalize() {};

ConcurrentClusterWorker::ConcurrentClusterWorker(
  ClusterAlgorithm *algo,
  std::istream &file,
  const int threads
) :
  ClusterWorker(file, threads),
  algo(algo) {};

void ConcurrentClusterWorker::operator()(size_t begin, size_t end) {
  DistanceElement d;
  std::vector<DistanceElement> buffer;
  buffer.reserve(100);
  std::cout << "Starting ConcurrentClusterWorker thread " << begin << std::endl;
  while (true) {
    mutex.lock();
    for (int i = 0; i < 100 && file; ++i) {
      file >> d;
      buffer.push_back(d);
    }
    mutex.unlock();
    if (buffer.size() == 0) break;
    for (auto d : buffer) {
      algo->operator()(d, begin);
    }
    buffer.clear();
  }

  mutex.lock();
  std::cout << "ConcurrentClusterWorker thread " << begin << " exiting" << std::endl;
  mutex.unlock();
}

void ConcurrentClusterWorker::finalize() {};

HierarchicalClusterWorker::HierarchicalClusterWorker(
  ClusterAlgorithm *algo,
  std::istream &file,
  const int threads,
  const int shards
) : ClusterWorker(file, threads), algo_list(shards), shards(shards),
thread_count(shards) {
  std::cout << "HieararchicalClusterWorker constructor start...";
  for (int i = 0; i < shards; ++i) {
    algo_list[i] = algo->make_child();
  }
  std::cout << "done" << std::endl;
}

void HierarchicalClusterWorker::operator()(size_t begin, size_t end) {
  DistanceElement d;
  std::vector<DistanceElement> buffer;
  buffer.reserve(100);
  size_t i = begin % shards;
  mutex.lock();
  ++thread_count[i];
  std::cout << "Starting HierarchicalClusterWorker thread " << begin
            << " (shard " << i
            << " with " << thread_count[i]
            << " previous threads)"
            << std::endl;
  mutex.unlock();
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
    mutex.lock();
    std::cout << "HierarchicalClusterWorker thread " << begin
              << " merging shard " << i
              << "..." << std::endl;
    mutex.unlock();
    algo_list[i]->merge_into_parent();
  }

    mutex.lock();
    std::cout << "HieararchicalClusterWorker thread " << begin
              << " exiting" << std::endl;
    mutex.unlock();
  }

  void HierarchicalClusterWorker::finalize() {

}
