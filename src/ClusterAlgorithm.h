#ifndef _CLUSTER_ALGORITHM_
#define _CLUSTER_ALGORITHM_

#include "single_linkage.h"
#include "DistanceConverter.h"
#include <vector>
#include <unordered_map>
#include <map>

class DistanceConsumer {
public:
  virtual void operator()(j_t seq1, j_t seq2, double dist)=0;
};

class ClusterAlgorithm : public DistanceConsumer {
protected:
  const DistanceConverter &dconv;
  const j_t n;
  const d_t m;
public:
  ClusterAlgorithm(const DistanceConverter &dconv, j_t n, d_t m) :
  dconv(dconv), n(n), m(m) {};

  // send consumer() pairwise distances to ensure it is up-to-date with this
  // clustering
  virtual void merge_into(DistanceConsumer &consumer) const=0;

  // send consumer() pairwise distances to ensure it is up-to-date with this
  // clustering
  virtual void merge_into(ClusterAlgorithm &consumer) const=0;

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
};

#endif
