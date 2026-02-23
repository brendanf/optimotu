#include "MappedClusterAlgorithm.h"

using MCA = MappedClusterAlgorithm;

std::vector<std::size_t> map_from_pair_generator(const PairGenerator & pg) {
  std::vector<std::size_t> fwd_map(pg.max_value());
  for (std::size_t i = 0; i < pg.max_value(); i++) {
    fwd_map[i] = pg.forward_map(i);
  }
  return fwd_map;
}

std::unordered_map<std::size_t, std::size_t> invert_map(const std::vector<std::size_t> &fwd_map) {
  std::unordered_map<std::size_t, std::size_t> rev_map;
  for (std::size_t i = 0; i < fwd_map.size(); i++) {
    rev_map[fwd_map[i]] = i;
  }
  return rev_map;
}

// Surrogate algorithm

MCA::Surrogate::Surrogate(
  ClusterAlgorithm * parent,
  const std::vector<std::size_t> &fwd_map
) :
  ClusterAlgorithm(parent),
  fwd_map(fwd_map)
{
}

void MCA::Surrogate::operator()(j_t seq1, j_t seq2, double dist, int thread) {
  (*parent)(fwd_map[seq1], fwd_map[seq2], dist, thread);
}
void MCA::Surrogate::operator()(j_t seq1, j_t seq2, int i, int thread) {
  (*parent)(fwd_map[seq1], fwd_map[seq2], i, thread);
}
void MCA::Surrogate::merge_into(DistanceConsumer &consumer) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::merge_into() is not implemented");
}
void MCA::Surrogate::merge_into(ClusterAlgorithm &consumer) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::merge_into() is not implemented");
}
void MCA::Surrogate::merge_into_parent() {
  OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::merge_into_parent() is not implemented");
}
ClusterAlgorithm * MCA::Surrogate::make_child() {
  OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::make_child() is not implemented");
}
ClusterAlgorithm * MCA::Surrogate::make_child(PairGenerator * pg) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::make_child() is not implemented");
}
double MCA::Surrogate::max_relevant(j_t seq1, j_t seq2, int thread) const {
  OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::max_relevant() is not implemented");
}
#ifdef OPTIMOTU_R
  Rcpp::List MCA::Surrogate::as_hclust(const Rcpp::CharacterVector &seqnames) const {
    OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::as_hclust() is not implemented");
  }
#endif // OPTIMOTU_R

// DistanceForwarder

MCA::DistanceForwarder::DistanceForwarder(
  const std::vector<std::size_t> &fwd_map,
  DistanceConsumer * consumer
) :
  fwd_map(fwd_map),
  target(consumer)
{}
void MCA::DistanceForwarder::operator()(PairGenerator & pg, double dist, int thread) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::DistanceForwarder::operator()() should not be called on a PairGenerator");
}
void MCA::DistanceForwarder::operator()(j_t seq1, j_t seq2, double dist, int thread) {
  (*target)(fwd_map[seq1], fwd_map[seq2], dist, thread);
}

// IndexForwarder

MCA::IndexForwarder::IndexForwarder(
  const std::vector<std::size_t> &fwd_map,
  ClusterAlgorithm * target
) :
  ClusterAlgorithm(target),
  fwd_map(fwd_map),
  target(target)
{}

// stub method
void MCA::IndexForwarder::operator()(PairGenerator & pg, double dist, int thread) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::operator()() should not be called on a PairGenerator");
}
void MCA::IndexForwarder::operator()(j_t seq1, j_t seq2, double dist, int thread) {
  (*target)(fwd_map[seq1], fwd_map[seq2], dist, thread);
}
void MCA::IndexForwarder::operator()(PairGenerator & pg, int i, int thread) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::operator()() should not be called on a PairGenerator");
}
void MCA::IndexForwarder::operator()(j_t seq1, j_t seq2, int i, int thread) {
  (*target)(fwd_map[seq1], fwd_map[seq2], i, thread);
}
void MCA::IndexForwarder::merge_into(DistanceConsumer &consumer) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::merge_into() is not implemented");
}
void MCA::IndexForwarder::merge_into(ClusterAlgorithm &consumer) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::merge_into() is not implemented");
}
void MCA::IndexForwarder::merge_into_parent() {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::merge_into_parent() is not implemented");
}

ClusterAlgorithm * MCA::IndexForwarder::make_child() {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::make_child() is not implemented");
}
ClusterAlgorithm * MCA::IndexForwarder::make_child(PairGenerator * pg) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::make_child() is not implemented");
}
double MCA::IndexForwarder::max_relevant(j_t seq1, j_t seq2, int thread) const {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::max_relevant() is not implemented");
}
#ifdef OPTIMOTU_R
  Rcpp::List MCA::IndexForwarder::as_hclust(const Rcpp::CharacterVector &seqnames) const {
    OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::as_hclust() is not implemented");
  }
#endif // OPTIMOTU_R


MCA::MappedClusterAlgorithm(SingleClusterAlgorithm * parent, PairGenerator * pg) :
  SingleClusterAlgorithm(parent, pg->max_value()),
  fwd_map(map_from_pair_generator(*pg)),
  rev_map(invert_map(fwd_map)),
  surrogate(std::make_unique<Surrogate>(parent, fwd_map)),
  inner(*(parent->make_inner_child(surrogate.get(), pg->max_value()))),
  pair_generator(pg)
{}

MCA::MappedClusterAlgorithm(
  SingleClusterAlgorithm * parent,
  SingleClusterAlgorithm * delegate,
  PairGenerator * pg
) :
  SingleClusterAlgorithm(parent, pg->max_value()),
  fwd_map(map_from_pair_generator(*pg)),
  rev_map(invert_map(fwd_map)),
  surrogate(std::make_unique<Surrogate>(delegate, fwd_map)),
  inner(*(parent->make_inner_child(surrogate.get(), pg->max_value()))),
  pair_generator(pg)
{}

MCA::MappedClusterAlgorithm(
  SingleClusterAlgorithm & inner,
  const std::vector<std::size_t> &fwd_map
) :
  SingleClusterAlgorithm(inner.dconv, inner.n),
  fwd_map(fwd_map),
  rev_map(invert_map(fwd_map)),
  inner(inner)
{}

MCA::MappedClusterAlgorithm(
  SingleClusterAlgorithm * parent,
  const std::vector<std::size_t> &fwd_map,
  const std::unordered_map<std::size_t, std::size_t> &rev_map
) :
  SingleClusterAlgorithm(parent->dconv, fwd_map.size()),
  fwd_map(fwd_map),
  rev_map(rev_map),
  surrogate(std::make_unique<Surrogate>(parent, fwd_map)),
  inner(*(parent->make_inner_child(surrogate.get(), fwd_map.size())))
{}

void MCA::operator()(j_t seq1, j_t seq2, double dist, int thread) {
  inner(rev_map.at(seq1), rev_map.at(seq2), dist, thread);
}
void MCA::operator()(j_t seq1, j_t seq2, int i, int thread) {
  inner(rev_map.at(seq1), rev_map.at(seq2), i, thread);
}

// If the PairGenerator is the same one used to create this object,
// which is the normal case, we can use the internal indices directly
// instead of round-tripping through the forward and reverse maps.
void MCA::operator()(PairGenerator & pg, double dist, int thread) {
  if (!pg) return;
  if (&pg == this->pair_generator) {
    inner(pg.i(), pg.j(), dist, thread);
  } else {
    inner(rev_map.at(pg.i0()), rev_map.at(pg.j0()), dist, thread);
  }
}
void MCA::operator()(PairGenerator & pg, int i, int thread) {
  if (!pg) return;
  if (&pg == this->pair_generator) {
    inner(pg.i(), pg.j(), i, thread);
  } else {
    inner(rev_map.at(pg.i0()), rev_map.at(pg.j0()), i, thread);
  }
}

void MCA::merge_into(DistanceConsumer &consumer) {
  DistanceForwarder distance_forwarder(fwd_map, &consumer);
  inner.merge_into(distance_forwarder);
}
void MCA::merge_into(ClusterAlgorithm &consumer) {
  IndexForwarder index_forwarder(fwd_map, &consumer);
  inner.merge_into(index_forwarder);
}

// This will actually cause the inner algorithm to merge its results into
// the surrogate, which maps the indices and forwards to the real parent.
void MCA::merge_into_parent() {
  inner.merge_into_parent();
}

// This is called by a SplitClusterWorker to create a child algorithm.
// If we are already a MappedClusterAlgorithm, that means we are part of a
// MultipleClusterAlgorithm, and we should return a new MappedClusterAlgorithm
// that wraps the inner algorithm's child.
SingleClusterAlgorithm * MCA::make_child() {
  auto child_ptr = new MappedClusterAlgorithm(&inner, fwd_map, rev_map);
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  children.push_back(std::move(child));
  return child_ptr;
}

/* This is called by a {Dist}SplitClusterWorker to create a child algorithm.
  If we are already a MappedClusterAlgorithm and we are now getting split,
  that means we are part of a MultipleClusterAlgorithm.

  This results in a "nested" MappedClusterAlgorithm, where a newly created
  "inner" MCA, which is a child of the current "outer" MCA, wraps a child of
  the current "outer" MCA's inner algorithm. However, calls do not pass through the "outer" MCA; the "inner" MCA should be called directly.

  The new "inner" MCA's reverse map maps from the global indices to the indices
  used by the "inner" MCA's inner algorithm. However, its forward map is *not*
  the inverse of the reverse map; instead, it maps from the indices used by the
  "inner" MCA's inner algorithm to its parent's indices, i.e. the indices used
  by the "outer" MCA's inner algorithm.  This is because the subset of
  sequences which defines the "outer" algorithm are being considered in order to optimize the clustering threshold(s) only within that subset, and we do not need to map its results back to global indices.
  */
MCA * MCA::make_child(PairGenerator * pg) {
  // Temporarily calculate a map from the innermost indices to the outermost
  // indices.
  std::vector<std::size_t> new_fwd_map;
  std::size_t i = 0, j = 0;
  while (i < pg->max_value() && j < fwd_map.size()) {
    if (pg->forward_map(i) == fwd_map[j]) {
      new_fwd_map.push_back(fwd_map[j]);
      ++i;
      ++j;
    } else if (pg->forward_map(i) < fwd_map[j]) {
      ++i;
    } else {
      ++j;
    }
  }
  // If there is no intersection between the subsets, return nullptr.
  if (new_fwd_map.size() == 0) {
    return nullptr;
  }

  auto new_rev_map = invert_map(new_fwd_map);
  // Now apply the old reverse map to the new forward map.
  for (size_t i = 0; i < new_fwd_map.size(); i++) {
    new_fwd_map[i] = rev_map.at(new_fwd_map[i]);
  }
  // This constructor causes our inner algorithm to create a child, which it
  // owns. Our child wraps the inner algorithm's child.
  auto child_ptr = new MappedClusterAlgorithm(&inner, new_fwd_map, new_rev_map);
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  children.push_back(std::move(child));
  return child_ptr;
}

// This should never be needed.
SingleClusterAlgorithm * MCA::make_inner_child(ClusterAlgorithm * parent, const j_t n) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::make_inner_child() is not implemented");
}

double MCA::max_relevant(j_t seq1, j_t seq2, int thread) const {
  return inner.max_relevant(rev_map.at(seq1), rev_map.at(seq2), thread);
}
double MCA::max_relevant(PairGenerator & pg, int thread) const {
  if (!pg) return -1.0;
  if (&pg == this->pair_generator) {
    return inner.max_relevant(pg.i(), pg.j(), thread);
  } else {
    return inner.max_relevant(rev_map.at(pg.i0()), rev_map.at(pg.j0()), thread);
  }
}

// If we are part of a MultipleClusterAlgorithm, then this makes sense.
// If we are part of a SplitClusterWorker, then this should not be called.
void MCA::write_to_matrix(internal_matrix_t &out) {
  inner.write_to_matrix(out);
}

#ifdef OPTIMOTU_R
  // This makes sense when we are part of a MultipleClusterAlgorithm,
  // but not when we are one part of a SplitClusterWorker.
  Rcpp::List MCA::as_hclust(const Rcpp::CharacterVector &seqnames) const {
    return inner.as_hclust(seqnames);
  }
#endif // OPTIMOTU_R
