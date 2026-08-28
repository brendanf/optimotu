// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_CLUSTERALGORITHM_H_INCLUDED
#define OPTIMOTU_CLUSTERALGORITHM_H_INCLUDED

#include "single_linkage.h"
#include "DistanceConverter.h"
#include "PairGenerator.h"
#include <vector>
#include <mutex>
#include <shared_mutex>
#include <iostream>
#include <deque>
#include <memory>
#include <algorithm>

#ifdef OPTIMOTU_R
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using internal_matrix_t = RcppParallel::RMatrix<int>;
using internal_matrix_ref_t = RcppParallel::RMatrix<int>;
using init_matrix_t = Rcpp::IntegerMatrix;
#else
using internal_matrix_t = std::vector<int>;
using internal_matrix_ref_t = std::vector<int>&;
using init_matrix_t = std::vector<int>&;
#endif //OPTIMOTU_R

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
  virtual void operator()(j_t seq1, j_t seq2, double dist, int thread = 0)=0;
  virtual void operator()(PairGenerator & pg, double dist, int thread = 0) {
    if (pg) {
      this->operator()(pg.i0(), pg.j0(), dist, thread);
    }
  }
  virtual void operator()(DistanceElement d, int thread = 0) {
    this->operator()(d.seq1, d.seq2, d.dist, thread);
  };
  virtual ~DistanceConsumer()=default;
};

class ClusterAlgorithm : public DistanceConsumer {
public:
  const DistanceConverter &dconv;
  const d_t m;
protected:
  mutable std::shared_timed_mutex mutex;
  ClusterAlgorithm * const parent = nullptr;
  std::deque<std::unique_ptr<ClusterAlgorithm>> children;

  // constructor for child objects
  ClusterAlgorithm(ClusterAlgorithm * parent) :
    dconv(parent->dconv), m(parent->m), parent(parent) {};

public:
  using DistanceConsumer::operator();

  virtual void finalize() {};

  // construct a ClusterAlgorithm with the given DistanceConverter
  ClusterAlgorithm(const DistanceConverter &dconv) :
  dconv(dconv),  m(dconv.m) {};

  // move constructor
  // ClusterAlgorithm(ClusterAlgorithm&& c) : dconv(c.dconv), m(c.m),
  // parent(c.parent), children(std::move(c.children)) {};

  // send consumer() pairwise distances to ensure it is up-to-date with this
  // clustering
  virtual void merge_into(DistanceConsumer &consumer)=0;

  // send consumer() pairwise distances to ensure it is up-to-date with this
  // clustering
  virtual void merge_into(ClusterAlgorithm &consumer)=0;

  // send parent pairwise distances to ensure it is up-to-date
  virtual void merge_into_parent() {
    if (parent) this->merge_into(*parent);
  }

  // Release a specific child algorithm previously created by make_child().
  // Safe to call when the child has finished merging into this parent.
  virtual void release_child(const ClusterAlgorithm * child_ptr) {
    std::unique_lock<std::shared_timed_mutex> lock(this->mutex, std::try_to_lock);
    if (!lock.owns_lock()) {
      return;
    }
    auto it = std::find_if(
      children.begin(),
      children.end(),
      [child_ptr](const std::unique_ptr<ClusterAlgorithm> &child) {
        return child.get() == child_ptr;
      }
    );
    if (it != children.end()) {
      children.erase(it);
    }
  }

  // create a copy of this algorithm, which will merge its results into this one
  // when it is finished
  virtual ClusterAlgorithm * make_child() = 0;
  virtual ClusterAlgorithm * make_child(PairGenerator * pg) = 0;
  // Child sized to fwd_map.size(); fwd_map[local] is the parent index.
  virtual ClusterAlgorithm *make_child(const std::vector<std::size_t> &fwd_map);

  // calculate the maximum distance between seq1 and seq2 which would actually
  // cause an update
  virtual double max_relevant(j_t seq1, j_t seq2, int thread = 0) const=0;
  virtual double max_relevant(PairGenerator & pg, int thread = 0) const {
    return this->max_relevant(pg.i0(), pg.j0(), thread);
  }

  // update clustering based on pairwise distance index between seq1 and seq2
  virtual void operator()(j_t seq1, j_t seq2, d_t i, int thread = 0)=0;
  virtual void operator()(PairGenerator & pg, d_t i, int thread = 0) {
    if (pg) {
      this->operator()(pg.i0(), pg.j0(), i, thread);
    }
  }

  // update clustering based on pairwise distance between seq1 and seq2
  void operator()(j_t seq1, j_t seq2, double dist, int thread = 0) override {
    if (seq1 == seq2) return;
    d_t i = dconv.convert(dist);
    (*this)(seq1, seq2, i);
  };

  virtual void operator()(PairGenerator & pg, double dist, int thread = 0) override {
    this->operator()(pg, dconv.convert(dist), thread);
  }

  // Whether this algorithm can safely consume pair updates in arbitrary order.
  // Algorithms that require strict within-pair or across-pairs ordering (such as
  // SLINK) should override and return false.
  virtual bool accepts_unordered_pairs() const {
    return true;
  }

#ifdef OPTIMOTU_R
  // convert the clustering results to an hclust object
  virtual Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const =0;
#endif //OPTIMOTU_R

};

template<int verbose> class MappedClusterAlgorithmImpl;
class SingleClusterAlgorithm : public ClusterAlgorithm {
  template<int verbose> friend class MappedClusterAlgorithmImpl;
public:
  const j_t n;
protected:

  // constructor for child objects
  SingleClusterAlgorithm(ClusterAlgorithm * parent, j_t n) :
    ClusterAlgorithm(parent), n(n) {};

  // Construct a SingleClusterAlgorithm with the given DistanceConverter, to
  // cluster n objects at m thresholds.
  SingleClusterAlgorithm(const DistanceConverter &dconv, const j_t n) :
    ClusterAlgorithm(dconv), n(n) {};

  // Construct a SingleClusterAlgorithm with the given DistanceConverter, to
  // cluster at m thresholds, with the number of objects determined by the size
  // of im. Optionally im can also be used as internal storage
  SingleClusterAlgorithm(const DistanceConverter &dconv, init_matrix_t im) :
    ClusterAlgorithm(dconv), n(im.size() / dconv.m) {};

  // move constructor
  // SingleClusterAlgorithm(SingleClusterAlgorithm&& c) : ClusterAlgorithm(c), n(c.n) {};

  // Called by MappedClusterAlgorithm to create the "inner" child algorithm
  // Should create a new SingleClusterAlgorithm with the same DistanceConverter
  // as this object, typically also of the same derived type. However, it
  // should assign the parent as requested. This object should own the new
  // child object!
  virtual SingleClusterAlgorithm * make_inner_child(
    ClusterAlgorithm * parent,
    const j_t n
  ) = 0;

public:

  // create a copy of this algorithm, which will merge its results into this one
  // when it is finished
  virtual SingleClusterAlgorithm * make_child() = 0;
  virtual SingleClusterAlgorithm * make_child(PairGenerator * pg) = 0;
  
  // write results to the provided clustering matrix
  virtual void write_to_matrix(internal_matrix_t &out)=0;
};

#endif //OPTIMOTU_CLUSTERALGORITHM_H_INCLUDED
