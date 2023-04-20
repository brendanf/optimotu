#ifndef _CLUSTER_ALGORITHM_
#define _CLUSTER_ALGORITHM_

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
  const j_t n;
  const d_t m;
  mutable tbb::queuing_rw_mutex mutex;
  ClusterAlgorithm * const parent = nullptr;
  bool own_child = false;
  std::set<ClusterAlgorithm *> children;

  ClusterAlgorithm(ClusterAlgorithm * parent) :
    dconv(parent->dconv), n(parent->n), m(parent->m), parent(parent) {};
public:
  using DistanceConsumer::operator();
  ClusterAlgorithm(const DistanceConverter &dconv, const j_t n, const d_t m) :
  dconv(dconv), n(n), m(m) {};

  ClusterAlgorithm(ClusterAlgorithm&& c) : dconv(c.dconv), n(c.n), m(c.m),
  parent(c.parent), own_child(c.own_child), children(std::move(c.children)){
    c.children.clear();
  };

  virtual ~ClusterAlgorithm() {
    for (auto c : children) {
      delete c;
    }
  };

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
  virtual void write_to_matrix(RcppParallel::RMatrix<int> &out)=0;
#endif

};

#endif
