#ifndef OPTIMOTU_CLUSTERALGORITHM_H_INCLUDED
#define OPTIMOTU_CLUSTERALGORITHM_H_INCLUDED

#include "single_linkage.h"
#include "DistanceConverter.h"
#include <vector>
#include <unordered_map>
#include <map>
#include <mutex>
#include <set>
#include <iostream>

#ifdef OPTIMOTU_R
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using internal_matrix_t = RcppParallel::RMatrix<int>;
using init_matrix_t = Rcpp::IntegerMatrix;
#else
using internal_matrix_t = std::vector<int>&;
using init_matrix_t = std::vector<int>&;
#endif

struct DistanceElement{
  std::size_t seq1, seq2;
  double dist;
  friend std::istream & operator >> (std::istream &in, DistanceElement &d);
};

inline std::istream & operator >> (std::istream &in, DistanceElement &d) {
  return in >> d.seq1 >> d.seq2 >> d.dist;
}

class DistanceConsumer {
public:
  virtual void operator()(j_t seq1, j_t seq2, double dist)=0;
  virtual void operator()(DistanceElement d) {
    this->operator()(d.seq1, d.seq2, d.dist);
  };
  virtual ~DistanceConsumer()=default;
};

class ClusterAlgorithm : public DistanceConsumer {
protected:
  const DistanceConverter &dconv;
  const d_t m;
  mutable tbb::queuing_rw_mutex mutex;
  ClusterAlgorithm * const parent = nullptr;
  bool own_child = false;
  std::deque<std::unique_ptr<ClusterAlgorithm>> children;

  // constructor for child objects
  ClusterAlgorithm(ClusterAlgorithm * parent) :
    dconv(parent->dconv), m(parent->m), parent(parent) {};
public:
  using DistanceConsumer::operator();

  // construct a ClusterAlgorithm with the given DistanceConverter
  ClusterAlgorithm(const DistanceConverter &dconv) :
  dconv(dconv),  m(dconv.m) {};

  // move constructor
  // ClusterAlgorithm(ClusterAlgorithm&& c) : dconv(c.dconv), m(c.m),
  // parent(c.parent), own_child(c.own_child), children(std::move(c.children)) {};

  // send consumer() pairwise distances to ensure it is up-to-date with this
  // clustering
  virtual void merge_into(DistanceConsumer &consumer)=0;

  // send consumer() pairwise distances to ensure it is up-to-date with this
  // clustering
  virtual void merge_into(ClusterAlgorithm &consumer)=0;

  // send parent pairwise distances to ensure it is up-to-date
  virtual void merge_into_parent() {
    if (parent && !own_child) this->merge_into(*parent);
  }

  // create a copy of this algorithm, which will merge its results into this one
  // when it is finished
  virtual ClusterAlgorithm * make_child() = 0;

  // calculate the maximum distance between seq1 and seq2 which would actually
  // cause an update
  virtual double max_relevant(j_t seq1, j_t seq2) const=0;

  // update clustering based on pairwise distance index between seq1 and seq2
  virtual void operator()(j_t seq1, j_t seq2, d_t i)=0;

  // update clustering based on pairwise distance between seq1 and seq2
  void operator()(j_t seq1, j_t seq2, double dist) override {
    if (seq1 == seq2) return;
    d_t i = dconv.convert(dist);
    (*this)(seq1, seq2, i);
  };

#ifdef OPTIMOTU_R
  // convert the clustering results to an hclust object
  virtual Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const =0;
#endif

};

class SingleClusterAlgorithm : public ClusterAlgorithm {
protected:
  const j_t n;

  // constructor for child objects
  SingleClusterAlgorithm(SingleClusterAlgorithm * parent) :
    ClusterAlgorithm(parent), n(parent->n) {};

  // construct a SingleClusterAlgorithm with the given DistanceConverter, to cluster
  // n objects at m thresholds.
  SingleClusterAlgorithm(const DistanceConverter &dconv, const j_t n) :
    ClusterAlgorithm(dconv), n(n) {};

  // construct a SingleClusterAlgorithm with the given DistanceConverter, to cluster
  // at m thresholds, with the number of objects determined by the size of im.
  // Optionally im can also be used as internal storage
  SingleClusterAlgorithm(const DistanceConverter &dconv, init_matrix_t im) :
    ClusterAlgorithm(dconv), n(im.size() / dconv.m) {};

  // move constructor
  // SingleClusterAlgorithm(SingleClusterAlgorithm&& c) : ClusterAlgorithm(c), n(c.n) {};

public:

  // write results to the provided clustering matrix
  virtual void write_to_matrix(internal_matrix_t &out)=0;

  // create a copy of this algorithm, which will merge its results into this one
  // when it is finished
  virtual SingleClusterAlgorithm * make_child() = 0;
};

#endif //OPTIMOTU_CLUSTERALGORITHM_H_INCLUDED
