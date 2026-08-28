#include <cassert>
#include <algorithm>
#include <limits>
#include <string>
#include <numeric>
#include "MultipleClusterAlgorithm.h"
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

// Base protected main constructor: initializer list only.
MCA::MultipleClusterAlgorithm(
  const ClusterAlgorithmFactory & factory,
  const std::vector<std::string> &names,
  const std::vector<std::vector<std::string>> &subset_names,
  const int threads,
  std::shared_ptr<MemoryBudgetTracker> memory_budget
) :
  ClusterAlgorithm(factory.dconv),
  factory(factory),
  names(names),
  subset_indices(calculate_subset_indices(names, subset_names)),
  subset_names(subset_names),
  threads(threads),
  memory_budget(memory_budget),
  subset_key(names.size()),
  fwd_map(subset_names.size()),
  subsets(),
  borrowed_subsets(),
  tracked_allocations()
{
  tracked_base_allocation = estimate_base_allocation(names, subset_names, threads);
  budget_acquire(this->memory_budget, tracked_base_allocation, "MultipleClusterAlgorithm base");
}

template<int verbose>
MultipleClusterAlgorithmImpl<verbose>::MultipleClusterAlgorithmImpl(MultipleClusterAlgorithm * parent) :
  MultipleClusterAlgorithm(parent) {}

template <int verbose>
MultipleClusterAlgorithmImpl<verbose>::MultipleClusterAlgorithmImpl(
    MultipleClusterAlgorithm *parent,
    PairGenerator *pg) : MultipleClusterAlgorithm(parent, pg) {}

template<int verbose>
MultipleClusterAlgorithmImpl<verbose>::MultipleClusterAlgorithmImpl(
  MultipleClusterAlgorithm * parent,
  const std::vector<std::string> names,
  const std::vector<std::vector<j_t>> subset_indices,
  const std::vector<std::vector<std::string>> subset_names,
  const std::vector<std::vector<j_t>> subset_key,
  const std::vector<std::unordered_map<j_t, j_t>> fwd_map,
  const std::vector<j_t> child_to_parent_map,
  PairGenerator * pg,
  const int threads,
  std::shared_ptr<MemoryBudgetTracker> memory_budget
) : MultipleClusterAlgorithm(parent, names, subset_indices, subset_names, subset_key, fwd_map, child_to_parent_map, pg, threads, memory_budget) {}

template<int verbose>
MultipleClusterAlgorithmImpl<verbose>::MultipleClusterAlgorithmImpl(
  const ClusterAlgorithmFactory & factory,
  const std::vector<std::string> &names,
  const std::vector<std::vector<std::string>> &subset_names,
  const int threads,
  int,
  std::shared_ptr<MemoryBudgetTracker> memory_budget
) :
  MultipleClusterAlgorithm(factory, names, subset_names, threads, memory_budget)
{
  OPTIMOTU_DEBUG(
    1,
    << "  Allocating " << subset_names.size() << " subsets..." << std::flush;
  );
  subsets.reserve(subset_names.size());

  OPTIMOTU_DEBUG(
    1,
    << "done\n  Generating sequence name key for " << names.size()
    << " sequence names..." << std::flush;
  );
  std::unordered_map<std::string, j_t> namekey;
  for (j_t i = 0; i < names.size(); i++) {
    namekey.emplace(names[i], i);
  }

  OPTIMOTU_DEBUG(
    1,
    << "done\n  Mapping " << subset_names.size()
    << " subsets..." << std::flush;
  );

  OPTIMOTU_DEBUG(
    2,
    << std::endl;
  );

  for (j_t i = 0; i < subset_names.size(); ++i) {
    const std::size_t subset_bytes = factory.estimate_bytes(subset_names[i].size());
    budget_acquire(this->memory_budget, subset_bytes, "MultipleClusterAlgorithm subset init");
    owned_subsets.emplace_back(factory.create(subset_names[i].size(), verbose));
    tracked_allocations.push_back(subset_bytes);
    subsets.push_back(owned_subsets.back().get());
    fwd_map[i].reserve(subset_names[i].size());
    OPTIMOTU_DEBUG(
      4,
      << "  Subset " << i << ":" << std::endl;
    );
    for (j_t j = 0; j < subset_names[i].size(); ++j) {
      auto f = namekey.find(subset_names[i][j]);
      assert(f != namekey.end());
      OPTIMOTU_DEBUG(
        4,
        << "    adding sequence " << f->second
        << " (" << subset_names[i][j]
        << ") to subset " << i
        << " at position " << fwd_map[i].size()
        << std::endl;
      );
      subset_key[f->second].push_back(i);
      fwd_map[i].emplace(f->second, j);
      OPTIMOTU_DEBUG(
        4,
        << "    sequence " << f->second
        << " now found in " << subset_key[j].size()
        << " subsets\n   precluster " << i
        << " now has " << fwd_map[i].size()
        << " sequences" << std::endl;
      );
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

MCA::MultipleClusterAlgorithm(MCA * parent) :
  ClusterAlgorithm(parent),
  factory(parent->factory),
  names(parent->names),
  subset_indices(parent->subset_indices),
  subset_names(parent->subset_names),
  threads(parent->threads),
  memory_budget(parent->memory_budget),
  subset_key(parent->subset_key),
  fwd_map(parent->fwd_map),
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
                                                                borrowed_subsets(),
                                                                tracked_allocations(),
                                                                tile_routing_parent(parent)
{
  subsets.resize(parent->subsets.size(), nullptr);

  std::unordered_set<std::size_t> tile_globals;
  tile_globals.reserve(pg->max_value() + 1);
  for (std::size_t i = 0; i <= pg->max_value(); ++i)
  {
    tile_globals.insert(pg->forward_map(i));
  }

  for (j_t j = 0; j < static_cast<j_t>(parent->subsets.size()); ++j)
  {
    std::vector<j_t> locals;
    locals.reserve(parent->fwd_map[j].size());
    for (const auto &kv : parent->fwd_map[j])
    {
      if (tile_globals.count(kv.first))
      {
        locals.push_back(kv.second);
      }
    }
    if (locals.size() < 2)
    {
      continue;
    }
    const std::size_t subset_bytes = factory.estimate_bytes(parent->subsets[j]->n);
    budget_acquire(
        memory_budget,
        subset_bytes,
        "MultipleClusterAlgorithm tile child subset");
    auto child_subset = parent->subsets[j]->make_child();
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

// This constructor is used to create a child algorithm for a
// MultipleClusterAlgorithm.
MCA::MultipleClusterAlgorithm(
  MCA * parent,
  const std::vector<std::string> names,
  const std::vector<std::vector<j_t>> subset_indices,
  const std::vector<std::vector<std::string>> subset_names,
  const std::vector<std::vector<j_t>> subset_key,
  const std::vector<std::unordered_map<j_t, j_t>> fwd_map,
  const std::vector<j_t> child_to_parent_map,
  PairGenerator * pg,
  const int threads,
  std::shared_ptr<MemoryBudgetTracker> memory_budget
) :
  ClusterAlgorithm{parent},
  factory{parent->factory},
  names{names},
  subset_indices{subset_indices},
  subset_names{subset_names},
  threads{threads},
  memory_budget{memory_budget ? memory_budget : parent->memory_budget},
  subset_key{subset_key},
  fwd_map{fwd_map},
  borrowed_subsets(),
  tracked_allocations()
  {
    tracked_base_allocation = estimate_base_allocation(names, subset_names, threads);
    budget_acquire(this->memory_budget, tracked_base_allocation, "MultipleClusterAlgorithm mapped child base");
    subsets.reserve(child_to_parent_map.size());
    for (const j_t j0 : child_to_parent_map) {
      const std::size_t subset_bytes = factory.estimate_bytes(subset_indices[subsets.size()].size());
      budget_acquire(this->memory_budget, subset_bytes, "MultipleClusterAlgorithm mapped child subset");
      auto child_subset = parent->subsets[j0]->make_child(pg);
      subsets.push_back(child_subset);
      borrowed_subsets.emplace_back(parent->subsets[j0], child_subset);
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
      filter_irrelevant &&
      !(dist < subsets[j]->max_relevant(s1, s2, thread))
    ) {
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
    double maxj = subsets[j]->max_relevant(
        it1->second,
        it2->second,
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

void MCA::write_to_matrix(std::vector<internal_matrix_t> &matrix_list) {
  for (size_t i = 0; i < this->subsets.size(); i++) {
    if (this->subsets[i])
    {
      this->subsets[i]->write_to_matrix(matrix_list[i]);
    }
  }
}

#ifdef OPTIMOTU_R
Rcpp::List MCA::as_hclust() const {
  std::vector<Rcpp::List> out;
  const auto &subset_name_list = routing_subset_names();
  for (size_t i = 0; i < this->subsets.size(); i++) {
    if (this->subsets[i])
    {
      out.push_back(
          this->subsets[i]->as_hclust(Rcpp::wrap(subset_name_list[i])));
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
