#include <cassert>
#include <algorithm>
#include <limits>
#include <string>
#include <numeric>
#include "MultipleClusterAlgorithm.h"
#include "MappedClusterAlgorithm.h"
#include "optimotu.h"

using MCA = MultipleClusterAlgorithm;

std::vector<std::vector<j_t>> calculate_subset_indices(
  const std::vector<std::string> &names,
  const std::vector<std::vector<std::string>> &subset_names
) {
  std::unordered_map<std::string, j_t> namekey;
  for (j_t i = 0; i < names.size(); i++) {
    namekey.emplace(names[i], i);
  }
  std::vector<std::vector<j_t>> out;
  out.reserve(subset_names.size());
  for (const auto & my_subset_names : subset_names) {
    std::vector<j_t> subset_indices;
    subset_indices.reserve(my_subset_names.size());
    for (const auto & name : my_subset_names) {
      auto f = namekey.find(name);
      if (f == namekey.end()) {
        OPTIMOTU_STOP("name '" + name +
           "' in subset is not found in the names of the sequences");
      }
      subset_indices.push_back(f->second);
    }
    std::sort(subset_indices.begin(), subset_indices.end());
    out.push_back(subset_indices);
  }
  return out;
}

namespace {

std::size_t estimate_base_allocation(
  const std::vector<std::string> &names,
  const std::vector<std::vector<std::string>> &subset_names,
  int threads
) {
  std::size_t total = sizeof(MultipleClusterAlgorithm);
  total += names.size() * sizeof(std::string);
  for (const auto &name : names) {
    total += name.size();
  }
  total += subset_names.size() * sizeof(std::vector<std::string>);
  for (const auto &subset : subset_names) {
    total += subset.size() * sizeof(std::string);
    for (const auto &name : subset) {
      total += name.size();
    }
  }
  total += names.size() * sizeof(std::vector<j_t>);
  total += subset_names.size() * sizeof(std::unordered_map<j_t, j_t>);
  total += static_cast<std::size_t>(threads) * sizeof(std::vector<j_t>);
  total += static_cast<std::size_t>(threads) * sizeof(std::pair<j_t, j_t>);
  return total;
}

void budget_acquire(
  const std::shared_ptr<MemoryBudgetTracker> &budget,
  std::size_t bytes,
  const std::string &context
) {
  if (budget) {
    budget->acquire(bytes, context);
  }
}

void budget_release(
  const std::shared_ptr<MemoryBudgetTracker> &budget,
  std::size_t bytes
) {
  if (budget) {
    budget->release(bytes);
  }
}

bool is_identity_perm(const std::vector<j_t> &perm)
{
  for (std::size_t i = 0; i < perm.size(); ++i)
  {
    if (perm[i] != static_cast<j_t>(i))
    {
      return false;
    }
  }
  return true;
}

#ifdef OPTIMOTU_R
void remap_hclust_tips(
    Rcpp::List hc,
    const std::vector<j_t> &sorted_to_which,
    const std::vector<std::string> &which_names)
{
  if (!is_identity_perm(sorted_to_which))
  {
    Rcpp::IntegerMatrix merge = hc["merge"];
    Rcpp::IntegerVector order = hc["order"];
    for (int i = 0; i < merge.nrow(); ++i)
    {
      for (int col = 0; col < 2; ++col)
      {
        int v = merge(i, col);
        if (v < 0)
        {
          const j_t sorted = static_cast<j_t>(-v - 1);
          merge(i, col) = -static_cast<int>(sorted_to_which[sorted]) - 1;
        }
      }
    }
    for (int i = 0; i < order.size(); ++i)
    {
      const j_t sorted = static_cast<j_t>(order[i] - 1);
      order[i] = static_cast<int>(sorted_to_which[sorted]) + 1;
    }
    hc["merge"] = merge;
    hc["order"] = order;
  }
  hc["labels"] = Rcpp::wrap(which_names);
}
#endif // OPTIMOTU_R

} // namespace

const std::vector<std::vector<j_t>> &MCA::routing_subset_key() const
{
  return tile_routing_parent ? tile_routing_parent->subset_key : subset_key;
}

const std::vector<std::unordered_map<j_t, j_t>> &MCA::routing_fwd_map() const
{
  return tile_routing_parent ? tile_routing_parent->fwd_map : fwd_map;
}

const std::vector<std::vector<std::string>> &MCA::routing_subset_names() const
{
  return tile_routing_parent ? tile_routing_parent->subset_names : subset_names;
}

const std::vector<std::string> &MCA::routing_names() const
{
  return tile_routing_parent ? tile_routing_parent->names : names;
}

const std::vector<std::vector<j_t>> &MCA::routing_subset_indices() const
{
  return tile_routing_parent ? tile_routing_parent->subset_indices
                             : subset_indices;
}

const std::vector<std::vector<j_t>> &MCA::routing_sorted_to_which() const
{
  return tile_routing_parent ? tile_routing_parent->sorted_to_which
                             : sorted_to_which;
}

// Base protected main constructor: initializer list only.
MCA::MultipleClusterAlgorithm(
    const ClusterAlgorithmFactory &factory,
    const std::vector<std::string> &names,
    const std::vector<std::vector<std::string>> &subset_names,
    const int threads,
    std::shared_ptr<MemoryBudgetTracker> memory_budget) : ClusterAlgorithm(factory.dconv),
                                                          factory(factory),
                                                          names(names),
                                                          subset_indices(calculate_subset_indices(names, subset_names)),
                                                          subset_names(subset_names),
                                                          threads(threads),
                                                          memory_budget(memory_budget),
                                                          subset_key(names.size()),
                                                          fwd_map(subset_names.size()),
                                                          sorted_to_which(subset_names.size()),
                                                          subsets(),
                                                          borrowed_subsets(),
                                                          tracked_allocations()
{
  tracked_base_allocation = estimate_base_allocation(names, subset_names, threads);
  budget_acquire(this->memory_budget, tracked_base_allocation, "MultipleClusterAlgorithm base");
}

template <int verbose>
MultipleClusterAlgorithmImpl<verbose>::MultipleClusterAlgorithmImpl(MultipleClusterAlgorithm *parent) : MultipleClusterAlgorithm(parent) {}

template <int verbose>
MultipleClusterAlgorithmImpl<verbose>::MultipleClusterAlgorithmImpl(
    MultipleClusterAlgorithm *parent,
    PairGenerator *pg) : MultipleClusterAlgorithm(parent, pg) {}

template <int verbose>
MultipleClusterAlgorithmImpl<verbose>::MultipleClusterAlgorithmImpl(
    const ClusterAlgorithmFactory &factory,
    const std::vector<std::string> &names,
    const std::vector<std::vector<std::string>> &subset_names,
    const int threads,
    int,
    std::shared_ptr<MemoryBudgetTracker> memory_budget) : MultipleClusterAlgorithm(factory, names, subset_names, threads, memory_budget)
{
  OPTIMOTU_DEBUG(
      1,
      << "  Allocating " << subset_names.size() << " subsets..." << std::flush;);
  subsets.reserve(subset_names.size());

  OPTIMOTU_DEBUG(
      1,
      << "done\n  Generating sequence name key for " << names.size()
      << " sequence names..." << std::flush;);
  std::unordered_map<std::string, j_t> namekey;
  for (j_t i = 0; i < names.size(); i++)
  {
    namekey.emplace(names[i], i);
  }

  OPTIMOTU_DEBUG(
      1,
      << "done\n  Mapping " << subset_names.size()
      << " subsets..." << std::flush;);

  OPTIMOTU_DEBUG(
      2,
      << std::endl;);

  for (j_t i = 0; i < subset_names.size(); ++i)
  {
    const std::size_t subset_bytes = factory.estimate_bytes(subset_names[i].size());
    budget_acquire(this->memory_budget, subset_bytes, "MultipleClusterAlgorithm subset init");
    owned_subsets.emplace_back(factory.create(subset_names[i].size(), verbose));
    tracked_allocations.push_back(subset_bytes);
    subsets.push_back(owned_subsets.back().get());
    fwd_map[i].reserve(subset_indices[i].size());
    sorted_to_which[i].assign(subset_names[i].size(), NO_CLUST);
    OPTIMOTU_DEBUG(
        4,
        << "  Subset " << i << ":" << std::endl;);
    for (j_t rank = 0; rank < subset_indices[i].size(); ++rank)
    {
      const j_t global = subset_indices[i][rank];
      OPTIMOTU_DEBUG(
          4,
          << "    adding sequence " << global
          << " (" << names[global]
          << ") to subset " << i
          << " at sorted position " << rank
          << std::endl;);
      subset_key[global].push_back(i);
      fwd_map[i].emplace(global, rank);
    }
    for (j_t which_pos = 0; which_pos < subset_names[i].size(); ++which_pos)
    {
      auto f = namekey.find(subset_names[i][which_pos]);
      assert(f != namekey.end());
      auto rank_it = fwd_map[i].find(f->second);
      assert(rank_it != fwd_map[i].end());
      sorted_to_which[i][rank_it->second] = which_pos;
    }
  }
}

// // Subset the names of the sequences to only the names of the sequences
// // which can be reached by the PairGenerator.
// std::vector<std::string> child_names(
//   const std::vector<std::string> & names,
//   const PairGenerator * const pg
// ) {
//   std::vector<std::string> out;
//   out.reserve(pg->max_value());
//   for (size_t i = 0; i < pg->max_value(); i++) {
//     out.push_back(names[pg->forward_map(i)]);
//   }
//   return out;
// }

// std::vector<std::vector<j_t>> child_subset_key(
//   const std::vector<std::vector<j_t>> & parent_subset_key,
//   const std::size_t parent_n_subsets,
//   const PairGenerator * const pg
// ) {
//   std::vector<std::vector<j_t>> out(parent_n_subsets);
//   for (std::size_t i = 0; i < pg->max_value(); i++) {
//     std::size_t i0 = pg->forward_map(i);
//   }
//   return out;
// }

// // Subset the indices of the sequences to only those which can be reached by
// // the PairGenerator.  This results in some empty subsets!
// std::vector<std::vector<j_t>> child_subset_indices(
//   const std::vector<std::vector<j_t>> & parent_subset_key,
//   const std::size_t parent_n_subsets,
//   const PairGenerator * const pg
// ) {
//   std::vector<std::vector<j_t>> out(parent_n_subsets);
//   // For each sequence that the PairGenerator can produce,
//   // find the subsets that it belongs to in the parent
//   // and add the mapped index to each of those subsets.
//   for (std::size_t i = 0; i < pg->max_value(); i++) {
//     std::size_t i0 = pg->forward_map(i);
//     for (const j_t j : parent_subset_key[i0]) {
//       out[j].push_back(i);
//     }
//   }
//   return out;
// }

// // Find the names of the sequences in each subset.
// std::vector<std::vector<std::string>> child_subset_names(
//   const std::vector<std::string> & child_names,
//   const std::vector<std::vector<j_t>> & child_subset_indices
// ) {
//   std::vector<std::vector<std::string>> out;
//   out.reserve(child_subset_indices.size());
//   for (const auto & my_subset : child_subset_indices) {
//     std::vector<std::string> my_names;
//     my_names.reserve(my_subset.size());
//     for (const auto & i : my_subset) {
//       my_names.push_back(child_names[i]);
//     }
//     out.push_back(my_names);
//   }
//   return out;
// }

MCA::MultipleClusterAlgorithm(MCA *parent) : ClusterAlgorithm(parent),
                                             factory(parent->factory),
                                             names(parent->names),
                                             subset_indices(parent->subset_indices),
                                             subset_names(parent->subset_names),
                                             threads(parent->threads),
                                             memory_budget(parent->memory_budget),
                                             subset_key(parent->subset_key),
                                             fwd_map(parent->fwd_map),
                                             sorted_to_which(parent->sorted_to_which),
                                             borrowed_subsets(),
                                             tracked_allocations()
{
  tracked_base_allocation = estimate_base_allocation(names, subset_names, threads);
  budget_acquire(this->memory_budget, tracked_base_allocation, "MultipleClusterAlgorithm child base");
  subsets.reserve(subset_names.size());
  for (auto & subset : parent->subsets) {
    const std::size_t subset_bytes = factory.estimate_bytes(subset->n);
    budget_acquire(this->memory_budget, subset_bytes, "MultipleClusterAlgorithm child subset");
    auto child_subset = subset->make_child();
    subsets.push_back(child_subset);
    borrowed_subsets.emplace_back(subset, child_subset);
    tracked_allocations.push_back(subset_bytes);
  }
}

MCA::MultipleClusterAlgorithm(MCA *parent, PairGenerator *pg) : ClusterAlgorithm(parent),
                                                                factory(parent->factory),
                                                                names(),
                                                                subset_indices(),
                                                                subset_names(),
                                                                threads(parent->threads),
                                                                memory_budget(parent->memory_budget),
                                                                subset_key(),
                                                                fwd_map(),
                                                                sorted_to_which(),
                                                                borrowed_subsets(),
                                                                tracked_allocations(),
                                                                tile_routing_parent(parent)
{
  subsets.resize(parent->subsets.size(), nullptr);
  tile_subset_locals.resize(parent->subsets.size());

  for (j_t j = 0; j < static_cast<j_t>(parent->subsets.size()); ++j)
  {
    auto intersection_fwd = subset_tile_fwd_map(parent->fwd_map[j], *pg);
    if (intersection_fwd.size() < 2)
    {
      continue;
    }
    for (std::size_t local : intersection_fwd)
    {
      tile_subset_locals[j].insert(static_cast<j_t>(local));
    }
    const std::size_t subset_bytes =
        factory.estimate_bytes(intersection_fwd.size());
    budget_acquire(
        memory_budget,
        subset_bytes,
        "MultipleClusterAlgorithm tile child subset");
    ClusterAlgorithm *child_subset = static_cast<ClusterAlgorithm *>(
                                         parent->subsets[j])
                                         ->make_child(intersection_fwd);
    if (!child_subset)
    {
      budget_release(memory_budget, subset_bytes);
      continue;
    }
    subsets[j] = static_cast<SingleClusterAlgorithm *>(child_subset);
    borrowed_subsets.emplace_back(parent->subsets[j], child_subset);
    tracked_allocations.push_back(subset_bytes);
  }
}

MCA::~MultipleClusterAlgorithm() {
  for (std::size_t bytes : tracked_allocations) {
    budget_release(memory_budget, bytes);
  }
  budget_release(memory_budget, tracked_base_allocation);
}

// Overlapping subsets for (seq1, seq2), cached per OS thread and MCA
// instance so max_relevant then apply_pair on the same pair skips a
// second set_intersection. Concurrent workers must not share this
// storage: they all pass thread=0, and tile indexes are not a valid
// thread id.
const std::vector<j_t> & MCA::ensure_whichsets(j_t seq1, j_t seq2) const {
  const auto &key = routing_subset_key();
  if (seq1 >= key.size() || seq2 >= key.size())
  {
    OPTIMOTU_STOP("MultipleClusterAlgorithm: sequence index out of range (seq1=" + std::to_string(seq1) + ", seq2=" + std::to_string(seq2) + ", names.size()=" + std::to_string(key.size()) + ")");
  }
  thread_local const MCA * cached_algo = nullptr;
  thread_local j_t cached_seq1 = std::numeric_limits<j_t>::max();
  thread_local j_t cached_seq2 = std::numeric_limits<j_t>::max();
  thread_local std::vector<j_t> whichsets;
  if (
    cached_algo == this &&
    cached_seq1 == seq1 &&
    cached_seq2 == seq2
  ) {
    return whichsets;
  }
  whichsets.clear();
  if (whichsets.capacity() < routing_subset_names().size())
  {
    whichsets.reserve(routing_subset_names().size());
  }
  std::set_intersection(
      key[seq1].begin(),
      key[seq1].end(),
      key[seq2].begin(),
      key[seq2].end(),
      std::back_inserter(whichsets));
  cached_algo = this;
  cached_seq1 = seq1;
  cached_seq2 = seq2;
  return whichsets;
}

// Forward a pair only to subsets for which it is still relevant. The worker
// threshold is max_relevant across all overlapping subsets, so without this
// filter a pair that is only needed by one subset is also applied to others
// whose own max_relevant is already below the distance. That changes
// single-linkage updates via convert/inverse rounding at cluster heights.
void MCA::apply_pair(
  j_t seq1,
  j_t seq2,
  d_t i,
  double dist,
  int thread,
  bool filter_irrelevant
) {
  if (seq1 == seq2) return;
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  const std::vector<j_t> & which = ensure_whichsets(seq1, seq2);
  const auto &maps = routing_fwd_map();
  for (j_t j : which) {
    if (!subsets[j])
    {
      continue;
    }
    const auto it1 = maps[j].find(seq1);
    const auto it2 = maps[j].find(seq2);
    if (it1 == maps[j].end() || it2 == maps[j].end())
    {
      continue;
    }
    j_t s1 = it1->second;
    j_t s2 = it2->second;
    if (
        tile_routing_parent &&
        subsets[j]->n < tile_routing_parent->subsets[j]->n)
    {
      if (
          j >= static_cast<j_t>(tile_subset_locals.size()) ||
          !tile_subset_locals[j].count(s1) ||
          !tile_subset_locals[j].count(s2))
      {
        continue;
      }
    }
    if (
        filter_irrelevant &&
        !(dist < subsets[j]->max_relevant(s1, s2, thread)))
    {
      continue;
    }
    (*subsets[j])(s1, s2, i, thread);
  }
}

void MCA::operator()(j_t seq1, j_t seq2, double dist, int thread) {
  // Sequence workers (including Edlib concurrent, which does not use a
  // PairGenerator) gate on max_relevant across overlapping subsets.
  apply_pair(seq1, seq2, dconv.convert(dist), dist, thread, true);
}

void MCA::operator()(j_t seq1, j_t seq2, int i, int thread) {
  apply_pair(seq1, seq2, i, dconv.inverse(i), thread, false);
}

void MCA::operator()(PairGenerator & pg, double dist, int thread) {
  if (!pg) return;
  // Sequence workers gate on max_relevant across all overlapping subsets.
  // Re-check each subset so a pair needed by one is not applied to others.
  apply_pair(
    pg.i0(),
    pg.j0(),
    dconv.convert(dist),
    dist,
    thread,
    true
  );
}

void MCA::operator()(DistanceElement d, int thread) {
  // Distance-matrix clustering sends every pair; ClusterTree decides what
  // is redundant. Independent distmx_cluster() of a subset also sees every
  // pair, so filtering here would make multi-subset results too sparse.
  apply_pair(d.seq1, d.seq2, dconv.convert(d.dist), d.dist, thread, false);
}

void MCA::finalize() {
  for (auto s : subsets)
  {
    if (s)
    {
      s->finalize();
    }
  }
}

// does not lock anything directly!
// relies on locks inside subset algorithms for thread safety.
double MCA::max_relevant(j_t seq1, j_t seq2, int thread) const {
  const auto &key = routing_subset_key();
  if (seq1 >= key.size() || seq2 >= key.size())
  {
    OPTIMOTU_STOP("MultipleClusterAlgorithm::max_relevant: sequence index out of range");
  }
  double max = -1.0;
  const std::vector<j_t> & which = ensure_whichsets(seq1, seq2);
  const auto &maps = routing_fwd_map();
  for (j_t j : which) {
    if (!subsets[j])
    {
      continue;
    }
    const auto it1 = maps[j].find(seq1);
    const auto it2 = maps[j].find(seq2);
    if (it1 == maps[j].end() || it2 == maps[j].end())
    {
      continue;
    }
    j_t s1 = it1->second;
    j_t s2 = it2->second;
    if (
        tile_routing_parent &&
        subsets[j]->n < tile_routing_parent->subsets[j]->n)
    {
      if (
          j >= static_cast<j_t>(tile_subset_locals.size()) ||
          !tile_subset_locals[j].count(s1) ||
          !tile_subset_locals[j].count(s2))
      {
        continue;
      }
    }
    double maxj = subsets[j]->max_relevant(
        s1,
        s2,
        thread);
    max = std::max(max, maxj);
  }
  return max;
}

void MCA::merge_into(DistanceConsumer &consumer) {
  for (auto & ss : this->subsets) {
    if (ss)
    {
      ss->merge_into(consumer);
    }
  }
}

void MCA::merge_into(ClusterAlgorithm &consumer) {
  if (!consumer.accepts_unordered_pairs()) {
    OPTIMOTU_STOP(
      "MultipleClusterAlgorithm::merge_into requires an unordered-pair consumer"
    );
  }
  for (auto & ss : this->subsets) {
    if (ss)
    {
      ss->merge_into(consumer);
    }
  }
}

// send consumer() pairwise distances to ensure it is up-to-date with this
// clustering
void MCA::merge_into(MCA &consumer) {
  for (size_t i = 0; i < this->subsets.size(); i++) {
    if (this->subsets[i])
    {
      this->subsets[i]->merge_into(*consumer.subsets[i]);
    }
  }
}

void MCA::merge_into_parent() {
  for (auto & ss : this->subsets) {
    if (ss)
    {
      ss->merge_into_parent();
    }
  }
}

MultipleClusterAlgorithm * MCA::make_child() {
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto child_ptr = new MultipleClusterAlgorithm(this);
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  this->children.push_back(std::move(child));
  return child_ptr;
}

template<int verbose>
MultipleClusterAlgorithm * MultipleClusterAlgorithmImpl<verbose>::make_child() {
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto child_ptr = new MultipleClusterAlgorithmImpl<verbose>(this);
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  this->children.push_back(std::move(child));
  return child_ptr;
}

MultipleClusterAlgorithm * MCA::make_child(PairGenerator * pg) {
  if (!pg)
  {
    return make_child();
  }
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto child_ptr = new MultipleClusterAlgorithm(this, pg);
  bool any_subset = false;
  for (const auto *subset : child_ptr->subsets)
  {
    if (subset)
    {
      any_subset = true;
      break;
    }
  }
  if (!any_subset)
  {
    delete child_ptr;
    return nullptr;
  }
  auto child = std::unique_ptr<ClusterAlgorithm>(
      (ClusterAlgorithm *)child_ptr);
  this->children.push_back(std::move(child));
  return child_ptr;
}

template<int verbose>
MultipleClusterAlgorithm * MultipleClusterAlgorithmImpl<verbose>::make_child(
  PairGenerator * pg
) {
  if (!pg)
  {
    return make_child();
  }
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto child_ptr = new MultipleClusterAlgorithmImpl<verbose>(this, pg);
  bool any_subset = false;
  for (const auto *subset : child_ptr->subsets)
  {
    if (subset)
    {
      any_subset = true;
      break;
    }
  }
  if (!any_subset)
  {
    delete child_ptr;
    return nullptr;
  }
  auto child = std::unique_ptr<ClusterAlgorithm>(
      (ClusterAlgorithm *)child_ptr);
  this->children.push_back(std::move(child));
  return child_ptr;
}

bool MCA::accepts_unordered_pairs() const {
  for (const auto * subset : subsets) {
    if (subset && !subset->accepts_unordered_pairs())
    {
      return false;
    }
  }
  return true;
}

void MCA::prepare_output()
{
  for (auto *s : this->subsets)
  {
    if (s)
    {
      s->prepare_output();
    }
  }
}

void MCA::write_threshold_row(std::size_t subset, d_t t, int *dest) const
{
  if (subset >= this->subsets.size() || !this->subsets[subset])
  {
    OPTIMOTU_STOP(
        "MultipleClusterAlgorithm::write_threshold_row: invalid subset");
  }
  const auto &perm = routing_sorted_to_which()[subset];
  const j_t nn = this->subsets[subset]->n;
  if (is_identity_perm(perm))
  {
    this->subsets[subset]->write_threshold_row(t, dest);
    return;
  }
  thread_local std::vector<int> tmp;
  tmp.resize(nn);
  this->subsets[subset]->write_threshold_row(t, tmp.data());
  for (j_t j = 0; j < nn; ++j)
  {
    dest[perm[j]] = tmp[j];
  }
}

void MCA::write_to_matrix(std::vector<internal_matrix_t> &matrix_list) {
  const auto &perm_all = routing_sorted_to_which();
  for (size_t i = 0; i < this->subsets.size(); i++) {
    if (!this->subsets[i])
    {
      continue;
    }
    const auto &perm = perm_all[i];
    if (is_identity_perm(perm))
    {
      this->subsets[i]->write_to_matrix(matrix_list[i]);
      continue;
    }
#ifdef OPTIMOTU_R
    Rcpp::IntegerMatrix tmp(matrix_list[i].nrow(), matrix_list[i].ncol());
    internal_matrix_t tmp_wrap(tmp);
    this->subsets[i]->write_to_matrix(tmp_wrap);
    const std::size_t m = matrix_list[i].nrow();
    const std::size_t n = matrix_list[i].ncol();
    for (std::size_t j = 0; j < n; ++j)
    {
      const std::size_t dest = perm[j];
      for (std::size_t r = 0; r < m; ++r)
      {
        matrix_list[i][dest * m + r] = tmp_wrap[j * m + r];
      }
    }
#else
    this->subsets[i]->write_to_matrix(matrix_list[i]);
#endif
  }
}

#ifdef OPTIMOTU_R
Rcpp::List MCA::as_hclust() const {
  std::vector<Rcpp::List> out;
  const auto &subset_name_list = routing_subset_names();
  const auto &all_names = routing_names();
  const auto &idx = routing_subset_indices();
  const auto &perm_all = routing_sorted_to_which();
  for (size_t i = 0; i < this->subsets.size(); i++) {
    if (this->subsets[i])
    {
      std::vector<std::string> sorted_labels;
      sorted_labels.reserve(idx[i].size());
      for (j_t global : idx[i])
      {
        sorted_labels.push_back(all_names[global]);
      }
      Rcpp::List hc = this->subsets[i]->as_hclust(Rcpp::wrap(sorted_labels));
      remap_hclust_tips(hc, perm_all[i], subset_name_list[i]);
      out.push_back(hc);
    }
  }
  return Rcpp::wrap(out);
}
Rcpp::List MCA::as_hclust(const Rcpp::CharacterVector &seqnames) const {
  return this->as_hclust();
}
#endif //OPTIMOTU_R

template class MultipleClusterAlgorithmImpl<0>;
template class MultipleClusterAlgorithmImpl<1>;
template class MultipleClusterAlgorithmImpl<2>;
template class MultipleClusterAlgorithmImpl<3>;
template class MultipleClusterAlgorithmImpl<4>;

#include <testthat.h>
#ifdef TESTTHAT_ENABLED
#include "ClusterAlgorithmFactory.h"

context("MCA subset routing")
{

  test_that("shuffled which maps to global ranks")
  {
    std::vector<std::string> names;
    for (int i = 0; i < 10; ++i)
    {
      names.push_back(std::to_string(i));
    }
    std::vector<std::vector<std::string>> subset{{"8", "0", "6", "2", "4"}};
    auto indices = calculate_subset_indices(names, subset);
    expect_true(indices.size() == 1u);
    expect_true(indices[0].size() == 5u);
    expect_true(indices[0][0] == 0);
    expect_true(indices[0][1] == 2);
    expect_true(indices[0][2] == 4);
    expect_true(indices[0][3] == 6);
    expect_true(indices[0][4] == 8);

    UniformDistanceConverter dconv(0.0, 0.4, 0.01);
    ClusterTreeFactory factory(dconv, 0);
    struct MCAPeek : public MultipleClusterAlgorithmImpl<0>
    {
      using MultipleClusterAlgorithmImpl<0>::MultipleClusterAlgorithmImpl;
      const std::vector<std::unordered_map<j_t, j_t>> &fwd() const
      {
        return fwd_map;
      }
      const std::vector<std::vector<j_t>> &perm() const
      {
        return sorted_to_which;
      }
    };
    MCAPeek algo(factory, names, subset, 1, 0);
    expect_true(algo.fwd().size() == 1u);
    expect_true(algo.fwd()[0].at(0) == 0);
    expect_true(algo.fwd()[0].at(2) == 1);
    expect_true(algo.fwd()[0].at(4) == 2);
    expect_true(algo.fwd()[0].at(6) == 3);
    expect_true(algo.fwd()[0].at(8) == 4);
    // which = {8, 0, 6, 2, 4}; sorted-local -> original which index
    expect_true(algo.perm()[0].size() == 5u);
    expect_true(algo.perm()[0][0] == 1);
    expect_true(algo.perm()[0][1] == 3);
    expect_true(algo.perm()[0][2] == 4);
    expect_true(algo.perm()[0][3] == 2);
    expect_true(algo.perm()[0][4] == 0);
  }
}
#endif
