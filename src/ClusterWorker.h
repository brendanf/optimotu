#ifndef OPTIMOTU_CLUSTERWORKER_H_INCLUDED
#define OPTIMOTU_CLUSTERWORKER_H_INCLUDED

#include <atomic>
#include <RcppParallel.h>
#include "ClusterAlgorithm.h"
#include "MultipleClusterAlgorithm.h"
#include "PairGenerator.h"
#include <unordered_map>
#include <memory>

class ClusterWorker : public RcppParallel::Worker {
protected:
  std::vector<std::unique_ptr<PairGenerator>> pair_generators;
  const int threads;
  int line_number = 0;
  ClusterAlgorithm * algo;
public:
  ClusterWorker(ClusterAlgorithm * algo, const int threads = 1);
  virtual ~ClusterWorker() = default;
  virtual void finalize() {};
  int n_threads();
};

template <typename distmx_t, typename id_t>
struct DistMatrixTranslator {
  DistMatrixTranslator(const distmx_t & distmx) {}
};

template <>
struct DistMatrixTranslator<Rcpp::DataFrame, int> {
  const Rcpp::IntegerVector _id1;
  const Rcpp::IntegerVector _id2;
  const Rcpp::NumericVector _dist;
  DistMatrixTranslator(const Rcpp::DataFrame & distmx) :
    _id1(distmx[0]), _id2(distmx[1]), _dist(distmx[2]) {}

  int id1(const int i) const {
    return _id1[i];
  }
  int id2(const int i) const {
    return _id2[i];
  }
  double dist(const int i) const {
    return _dist[i];
  }
};

template <>
struct DistMatrixTranslator<Rcpp::DataFrame, std::string> {
  const Rcpp::CharacterVector _id1;
  const Rcpp::CharacterVector _id2;
  const Rcpp::NumericVector _dist;
  DistMatrixTranslator(const Rcpp::DataFrame & distmx) :
    _id1(distmx[0]), _id2(distmx[1]), _dist(distmx[2]) {}
  std::string id1(const int i) const {
    Rcpp::String s = _id1[i];
    return std::string(s.get_cstring());
  }
  std::string id2(const int i) const {
    Rcpp::String s = _id2[i];
    return std::string(s.get_cstring());
  }
  double dist(const int i) const {
    return _dist[i];
  }
};

template <typename distmx_t, typename id_t>
class ClusterWorkerImpl : public ClusterWorker {
protected:
  distmx_t &distmx;
  typedef std::unordered_map<id_t, int> id_map_type;
  std::shared_ptr<id_map_type> id_map;
  std::mutex mutex;

  DistMatrixTranslator<distmx_t, id_t> translator;
  // sets d to the next line of the distance matrix
  // returns false if there are no more lines
  bool next_line(DistanceElement & d);
public:
  typedef std::vector<id_t> id_list_type;
  ClusterWorkerImpl(
    ClusterAlgorithm * algo,
    distmx_t &distmx,
    const int threads = 1
  );
  ClusterWorkerImpl(
    ClusterAlgorithm * algo,
    distmx_t &distmx,
    const id_list_type & id_list,
    const int threads = 1
  );
};

// The MergeClusterWorker is a ClusterWorker that works on separate
// instances of the ClusterAlgorithm, one for each thread.
// This means that the individual instances can be smaller, and there
// is no worry about concurrency issues among the different instances.
// They do each need to merge their results into the parent algorithm
// when they are finished, but this can typically be done using
// fewer steps than the original clustering.
// Implementations need to use PairGenerator::i() and PairGenerator::j()
// to get the indices of the pair in the cluster algorithm, but
// PairGenerator::i0() and PairGenerator::j0() to get the indices of the pair
// in the original sequence list.
template <class distmx_t, class id_t>
class MergeClusterWorker : public ClusterWorkerImpl<distmx_t, id_t> {
protected:
  using ClusterWorker::algo;
  using ClusterWorker::threads;
  std::vector<ClusterAlgorithm*> algo_list;
  using ClusterWorkerImpl<distmx_t, id_t>::next_line;
public:
  using typename ClusterWorkerImpl<distmx_t, id_t>::id_list_type;
  MergeClusterWorker(
    ClusterAlgorithm * algo,
    distmx_t &distmx,
    const int threads = 1
    );
  MergeClusterWorker(
    ClusterAlgorithm * algo,
    distmx_t &distmx,
    const id_list_type & id_list,
    const int threads = 1
  );

  void operator()(size_t begin, size_t end) override;
};

// The ConcurrentClusterWorker is a ClusterWorker that works on a single
// instance of the ClusterAlgorithm, but uses multiple threads to process
// the pairs.
// Implementations need to use PairGenerator::i0() and PairGenerator::j0()
// for indexes in both the cluster algorithm and the original sequence list.
template <class distmx_t, class id_t>
class ConcurrentClusterWorker : public ClusterWorkerImpl<distmx_t, id_t> {
protected:
  using ClusterWorker::algo;
  using ClusterWorker::threads;
  using ClusterWorkerImpl<distmx_t, id_t>::next_line;
public:
  using typename ClusterWorkerImpl<distmx_t, id_t>::id_list_type;
  using ClusterWorkerImpl<distmx_t, id_t>::ClusterWorkerImpl;

  void operator()(size_t begin, size_t end) override;
};

// The HierarchicalClusterWorker is a ClusterWorker that works on
// multiple instances of the ClusterAlgorithm (shards), each working on
// multiple subsets of pairs (threads). This means that the different instances
// can be smaller than in the serial case or the ConcurrentClusterWorker,
// but there is still some capability for concurrency.
// Unfortunately there are potentially three different indexes to keep track of:
// - the index of the pair in the pair generator (thread-level) -- i() and j()
// - the index of the pair in the cluster algorithm (shard-level) -- ??
// - the index of the pair in the original sequence list (global) -- i0() and j0()
// The first and last are easy to keep track of, but the middle is not implemented.
template <class distmx_t, class id_t>
class HierarchicalClusterWorker : public ClusterWorkerImpl<distmx_t, id_t> {
protected:
  using ClusterWorker::algo;
  using ClusterWorker::threads;
  std::vector<ClusterAlgorithm*> algo_list;
  std::vector<std::atomic_size_t> thread_count;
  const int shards;
  using ClusterWorkerImpl<distmx_t, id_t>::next_line;
public:
  using typename ClusterWorkerImpl<distmx_t, id_t>::id_list_type;

  HierarchicalClusterWorker(
    ClusterAlgorithm *algo,
    distmx_t & distmx,
    const int threads,
    const int shards
  );

  HierarchicalClusterWorker(
    ClusterAlgorithm *algo,
    distmx_t & distmx,
    const id_list_type & id_list,
    const int threads,
    const int shards
  );

  void operator()(size_t begin, size_t end) override;
};

#endif //OPTIMOTU_CLUSTERWORKER_H_INCLUDED
