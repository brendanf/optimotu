#include "MappedClusterAlgorithm.h"
#include "optimotu.h"

ClusterAlgorithm *ClusterAlgorithm::make_child(
    const std::vector<std::size_t> &fwd_map)
{
  (void)fwd_map;
  OPTIMOTU_STOP("ClusterAlgorithm::make_child(fwd_map) is not implemented");
  return nullptr;
}

using MCA = MappedClusterAlgorithm;

std::vector<std::size_t> map_from_pair_generator(const PairGenerator & pg) {
  std::vector<std::size_t> fwd_map(pg.max_value() + 1);
  for (std::size_t i = 0; i <= pg.max_value(); i++) {
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

// ---------------------------------------------------------------------------
// MappedClusterAlgorithm base
// ---------------------------------------------------------------------------

MCA::MappedClusterAlgorithm(ClusterAlgorithm * parent, j_t n) :
  SingleClusterAlgorithm(parent, n)
{}

MCA::MappedClusterAlgorithm(const DistanceConverter & dconv, j_t n) :
  SingleClusterAlgorithm(dconv, n)
{}

void MCA::merge_into_parent() {
  get_inner().merge_into_parent();
}

void MCA::write_to_matrix(internal_matrix_t &out) {
  get_inner().write_to_matrix(out);
}

bool MCA::accepts_unordered_pairs() const {
  return get_inner().accepts_unordered_pairs();
}

#ifdef OPTIMOTU_R
Rcpp::List MCA::as_hclust(const Rcpp::CharacterVector &seqnames) const {
  return get_inner().as_hclust(seqnames);
}
#endif // OPTIMOTU_R

// ---------------------------------------------------------------------------
// MappedClusterAlgorithmImpl<verbose>::Surrogate
// ---------------------------------------------------------------------------

template<int verbose>
MappedClusterAlgorithmImpl<verbose>::Surrogate::Surrogate(
  ClusterAlgorithm* parent,
  const std::vector<std::size_t>& fwd_map
) :
  ClusterAlgorithm(parent),
  fwd_map(fwd_map)
{}

template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::Surrogate::operator()(j_t seq1, j_t seq2, double dist, int thread) {
  if (
      seq1 < 0 || seq2 < 0 ||
      static_cast<std::size_t>(seq1) >= fwd_map.size() ||
      static_cast<std::size_t>(seq2) >= fwd_map.size())
  {
    return;
  }
  (*parent)(fwd_map[seq1], fwd_map[seq2], dist, thread);
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::Surrogate::operator()(j_t seq1, j_t seq2, int i, int thread) {
  if (
      seq1 < 0 || seq2 < 0 ||
      static_cast<std::size_t>(seq1) >= fwd_map.size() ||
      static_cast<std::size_t>(seq2) >= fwd_map.size())
  {
    return;
  }
  (*parent)(fwd_map[seq1], fwd_map[seq2], i, thread);
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::Surrogate::merge_into(DistanceConsumer &consumer) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::merge_into() is not implemented");
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::Surrogate::merge_into(ClusterAlgorithm &consumer) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::merge_into() is not implemented");
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::Surrogate::merge_into_parent() {
  OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::merge_into_parent() is not implemented");
}
template<int verbose>
ClusterAlgorithm* MappedClusterAlgorithmImpl<verbose>::Surrogate::make_child() {
  OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::make_child() is not implemented");
}
template<int verbose>
ClusterAlgorithm* MappedClusterAlgorithmImpl<verbose>::Surrogate::make_child(PairGenerator* pg) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::make_child() is not implemented");
}
template<int verbose>
double MappedClusterAlgorithmImpl<verbose>::Surrogate::max_relevant(j_t seq1, j_t seq2, int thread) const {
  OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::max_relevant() is not implemented");
}
template<int verbose>
bool MappedClusterAlgorithmImpl<verbose>::Surrogate::accepts_unordered_pairs() const {
  return parent->accepts_unordered_pairs();
}
#ifdef OPTIMOTU_R
template<int verbose>
Rcpp::List MappedClusterAlgorithmImpl<verbose>::Surrogate::as_hclust(const Rcpp::CharacterVector &seqnames) const {
  OPTIMOTU_STOP("MappedClusterAlgorithm::Surrogate::as_hclust() is not implemented");
}
#endif

// ---------------------------------------------------------------------------
// MappedClusterAlgorithmImpl<verbose>::DistanceForwarder
// ---------------------------------------------------------------------------

template<int verbose>
MappedClusterAlgorithmImpl<verbose>::DistanceForwarder::DistanceForwarder(
  const std::vector<std::size_t>& fwd_map,
  DistanceConsumer* consumer
) :
  fwd_map(fwd_map),
  target(consumer)
{}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::DistanceForwarder::operator()(PairGenerator& pg, double dist, int thread) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::DistanceForwarder::operator()() should not be called on a PairGenerator");
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::DistanceForwarder::operator()(j_t seq1, j_t seq2, double dist, int thread) {
  if (
      seq1 < 0 || seq2 < 0 ||
      static_cast<std::size_t>(seq1) >= fwd_map.size() ||
      static_cast<std::size_t>(seq2) >= fwd_map.size())
  {
    return;
  }
  (*target)(fwd_map[seq1], fwd_map[seq2], dist, thread);
}

// ---------------------------------------------------------------------------
// MappedClusterAlgorithmImpl<verbose>::IndexForwarder
// ---------------------------------------------------------------------------

template<int verbose>
MappedClusterAlgorithmImpl<verbose>::IndexForwarder::IndexForwarder(
  const std::vector<std::size_t>& fwd_map,
  ClusterAlgorithm* target
) :
  ClusterAlgorithm(target),
  fwd_map(fwd_map),
  target(target)
{}

template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::IndexForwarder::operator()(PairGenerator& pg, double dist, int thread) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::operator()() should not be called on a PairGenerator");
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::IndexForwarder::operator()(j_t seq1, j_t seq2, double dist, int thread) {
  if (
      seq1 < 0 || seq2 < 0 ||
      static_cast<std::size_t>(seq1) >= fwd_map.size() ||
      static_cast<std::size_t>(seq2) >= fwd_map.size())
  {
    return;
  }
  (*target)(fwd_map[seq1], fwd_map[seq2], dist, thread);
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::IndexForwarder::operator()(PairGenerator& pg, int i, int thread) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::operator()() should not be called on a PairGenerator");
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::IndexForwarder::operator()(j_t seq1, j_t seq2, int i, int thread) {
  if (
      seq1 < 0 || seq2 < 0 ||
      static_cast<std::size_t>(seq1) >= fwd_map.size() ||
      static_cast<std::size_t>(seq2) >= fwd_map.size())
  {
    return;
  }
  (*target)(fwd_map[seq1], fwd_map[seq2], i, thread);
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::IndexForwarder::merge_into(DistanceConsumer& consumer) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::merge_into() is not implemented");
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::IndexForwarder::merge_into(ClusterAlgorithm& consumer) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::merge_into() is not implemented");
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::IndexForwarder::merge_into_parent() {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::merge_into_parent() is not implemented");
}
template<int verbose>
ClusterAlgorithm* MappedClusterAlgorithmImpl<verbose>::IndexForwarder::make_child() {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::make_child() is not implemented");
}
template<int verbose>
ClusterAlgorithm* MappedClusterAlgorithmImpl<verbose>::IndexForwarder::make_child(PairGenerator* pg) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::make_child() is not implemented");
}
template<int verbose>
double MappedClusterAlgorithmImpl<verbose>::IndexForwarder::max_relevant(j_t seq1, j_t seq2, int thread) const {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::max_relevant() is not implemented");
}
template<int verbose>
bool MappedClusterAlgorithmImpl<verbose>::IndexForwarder::accepts_unordered_pairs() const {
  return target->accepts_unordered_pairs();
}
#ifdef OPTIMOTU_R
template<int verbose>
Rcpp::List MappedClusterAlgorithmImpl<verbose>::IndexForwarder::as_hclust(const Rcpp::CharacterVector& seqnames) const {
  OPTIMOTU_STOP("MappedClusterAlgorithm::IndexForwarder::as_hclust() is not implemented");
}
#endif

// ---------------------------------------------------------------------------
// MappedClusterAlgorithmImpl<verbose>::make_inner_child
// ---------------------------------------------------------------------------

template<int verbose>
SingleClusterAlgorithm* MappedClusterAlgorithmImpl<verbose>::make_inner_child(ClusterAlgorithm* parent, const j_t n) {
  OPTIMOTU_STOP("MappedClusterAlgorithm::make_inner_child() is not implemented");
}

// ---------------------------------------------------------------------------
// MappedClusterAlgorithmImpl<verbose> constructors
// ---------------------------------------------------------------------------

template<int verbose>
MappedClusterAlgorithmImpl<verbose>::MappedClusterAlgorithmImpl(SingleClusterAlgorithm* parent, PairGenerator* pg) :
  MappedClusterAlgorithm(parent, pg->max_value() + 1),
  fwd_map(map_from_pair_generator(*pg)),
  rev_map(invert_map(fwd_map)),
  surrogate(std::make_unique<Surrogate>(parent, fwd_map)),
  inner(parent->make_inner_child(surrogate.get(), pg->max_value() + 1)),
  pair_generator(pg)
{}

template<int verbose>
MappedClusterAlgorithmImpl<verbose>::MappedClusterAlgorithmImpl(
  SingleClusterAlgorithm* parent,
  SingleClusterAlgorithm* delegate,
  PairGenerator* pg
) :
  MappedClusterAlgorithm(parent, pg->max_value() + 1),
  fwd_map(map_from_pair_generator(*pg)),
  rev_map(invert_map(fwd_map)),
  surrogate(std::make_unique<Surrogate>(delegate, fwd_map)),
  inner(parent->make_inner_child(surrogate.get(), pg->max_value() + 1)),
  pair_generator(pg)
{}

template <int verbose>
MappedClusterAlgorithmImpl<verbose>::MappedClusterAlgorithmImpl(
    SingleClusterAlgorithm *parent,
    SingleClusterAlgorithm *delegate,
    const std::vector<std::size_t> &fwd_map,
    const std::unordered_map<std::size_t, std::size_t> &rev_map) : MappedClusterAlgorithm(parent, fwd_map.size()),
                                                                   fwd_map(fwd_map),
                                                                   rev_map(rev_map),
                                                                   surrogate(std::make_unique<Surrogate>(delegate, fwd_map)),
                                                                   inner(parent->make_inner_child(surrogate.get(), fwd_map.size())),
                                                                   pair_generator(nullptr)
{
}

template<int verbose>
MappedClusterAlgorithmImpl<verbose>::MappedClusterAlgorithmImpl(
  SingleClusterAlgorithm& inner_ref, const std::vector<std::size_t>& fwd_map
) :
  MappedClusterAlgorithm(inner_ref.dconv, inner_ref.n),
  fwd_map(fwd_map),
  rev_map(invert_map(fwd_map)),
  surrogate(),
  inner(&inner_ref),
  pair_generator(nullptr)
{}

template<int verbose>
MappedClusterAlgorithmImpl<verbose>::MappedClusterAlgorithmImpl(
  SingleClusterAlgorithm* parent,
  const std::vector<std::size_t>& fwd_map,
  const std::unordered_map<std::size_t, std::size_t>& rev_map
) :
  MappedClusterAlgorithm(parent->dconv, fwd_map.size()),
  fwd_map(fwd_map),
  rev_map(rev_map),
  surrogate(std::make_unique<Surrogate>(parent, fwd_map)),
  inner(parent->make_inner_child(surrogate.get(), fwd_map.size())),
  pair_generator(nullptr)
{}

// ---------------------------------------------------------------------------
// MappedClusterAlgorithmImpl<verbose>::operator(), merge_into, make_child
// ---------------------------------------------------------------------------

template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::operator()(j_t seq1, j_t seq2, double dist, int thread) {
  (*inner)(rev_map.at(seq1), rev_map.at(seq2), dist, thread);
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::operator()(j_t seq1, j_t seq2, int i, int thread) {
  (*inner)(rev_map.at(seq1), rev_map.at(seq2), i, thread);
}

template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::operator()(PairGenerator& pg, double dist, int thread) {
  if (!pg) return;
  if (&pg == pair_generator) {
    (*inner)(pg.i(), pg.j(), dist, thread);
  } else {
    (*inner)(rev_map.at(pg.i0()), rev_map.at(pg.j0()), dist, thread);
  }
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::operator()(PairGenerator& pg, int i, int thread) {
  if (!pg) return;
  if (&pg == pair_generator) {
    (*inner)(pg.i(), pg.j(), i, thread);
  } else {
    (*inner)(rev_map.at(pg.i0()), rev_map.at(pg.j0()), i, thread);
  }
}

template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::merge_into(DistanceConsumer& consumer) {
  DistanceForwarder distance_forwarder(fwd_map, &consumer);
  inner->merge_into(distance_forwarder);
}
template<int verbose>
void MappedClusterAlgorithmImpl<verbose>::merge_into(ClusterAlgorithm& consumer) {
  if (inner->accepts_unordered_pairs() && !consumer.accepts_unordered_pairs()) {
    OPTIMOTU_STOP(
      "MappedClusterAlgorithm::merge_into requires an unordered-pair consumer "
      "for unordered inner algorithms"
    );
  }
  IndexForwarder index_forwarder(fwd_map, &consumer);
  inner->merge_into(index_forwarder);
}

template<int verbose>
SingleClusterAlgorithm* MappedClusterAlgorithmImpl<verbose>::make_child() {
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto child_ptr = new MappedClusterAlgorithmImpl<verbose>(*inner, fwd_map);
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  this->children.push_back(std::move(child));
  return child_ptr;
}

template<int verbose>
MappedClusterAlgorithm* MappedClusterAlgorithmImpl<verbose>::make_child(PairGenerator* pg) {
  std::vector<std::size_t> new_fwd_map;
  std::size_t i = 0, j = 0;
  while (i <= pg->max_value() && j < fwd_map.size()) {
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
  if (new_fwd_map.size() == 0) {
    return nullptr;
  }
  auto new_rev_map = invert_map(new_fwd_map);
  for (size_t k = 0; k < new_fwd_map.size(); k++) {
    new_fwd_map[k] = rev_map.at(new_fwd_map[k]);
  }
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto child_ptr = new MappedClusterAlgorithmImpl<verbose>(*inner, new_fwd_map);
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  this->children.push_back(std::move(child));
  return child_ptr;
}

template<int verbose>
double MappedClusterAlgorithmImpl<verbose>::max_relevant(j_t seq1, j_t seq2, int thread) const {
  return inner->max_relevant(rev_map.at(seq1), rev_map.at(seq2), thread);
}

template<int verbose>
double MappedClusterAlgorithmImpl<verbose>::max_relevant(PairGenerator& pg, int thread) const {
  if (!pg) return -1.0;
  if (&pg == pair_generator) {
    return inner->max_relevant(pg.i(), pg.j(), thread);
  } else {
    return inner->max_relevant(rev_map.at(pg.i0()), rev_map.at(pg.j0()), thread);
  }
}

template class MappedClusterAlgorithmImpl<0>;
template class MappedClusterAlgorithmImpl<1>;
template class MappedClusterAlgorithmImpl<2>;
template class MappedClusterAlgorithmImpl<3>;
template class MappedClusterAlgorithmImpl<4>;