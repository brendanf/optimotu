#ifndef OPTIMOTU_CLUSTERWORKER_H_INCLUDED
#define OPTIMOTU_CLUSTERWORKER_H_INCLUDED

#include <atomic>
#include <RcppParallel.h>
#include "ClusterAlgorithm.h"
#include "MultipleClusterAlgorithm.h"
#include <unordered_map>
#include <memory>

class ClusterWorker : public RcppParallel::Worker {
protected:
  const int threads;
  int line_number = 0;
  ClusterAlgorithm * algo;
public:
  ClusterWorker(ClusterAlgorithm * algo, const int threads);
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
  ClusterWorkerImpl(ClusterAlgorithm * algo, distmx_t &distmx, const int threads);
  ClusterWorkerImpl(ClusterAlgorithm * algo, distmx_t &distmx,
                    const id_list_type & id_list, const int threads);
};

template <class distmx_t, class id_t>
class MergeClusterWorker : public ClusterWorkerImpl<distmx_t, id_t> {
protected:
  using ClusterWorker::algo;
  using ClusterWorker::threads;
  std::vector<ClusterAlgorithm*> algo_list;
  using ClusterWorkerImpl<distmx_t, id_t>::next_line;
public:
  using typename ClusterWorkerImpl<distmx_t, id_t>::id_list_type;
  MergeClusterWorker(ClusterAlgorithm * algo, distmx_t &distmx, const int threads);
  MergeClusterWorker(ClusterAlgorithm * algo, distmx_t &distmx,
                    const id_list_type & id_list, const int threads);

  void operator()(size_t begin, size_t end) override;
};

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
