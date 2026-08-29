// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "ClusterSLINK.h"
#include <algorithm>
#include <cstdint>

template <int verbose>
ClusterSLINK<verbose>::ClusterSLINK(ClusterAlgorithm *parent, j_t n) : SingleClusterAlgorithm(parent, n) {}

template <int verbose>
ClusterSLINK<verbose>::ClusterSLINK(const DistanceConverter &dconv, const j_t n) : SingleClusterAlgorithm(dconv, n) {}

template <int verbose>
ClusterSLINK<verbose>::ClusterSLINK(const DistanceConverter &dconv, init_matrix_t im) : SingleClusterAlgorithm(dconv, im) {}

template <int verbose>
void ClusterSLINK<verbose>::ensure_slink()
{
  if (!slink)
  {
    slink = std::make_unique<SlinkState>();
    slink->Pi.resize(n);
    slink->Lambda.assign(n, m);
    slink->M.assign(n, m);
  }
}

template <int verbose>
void ClusterSLINK<verbose>::init_iter()
{
  for (j_t k = 0; k < slink->slink_seq1; k++)
  {
    slink->M[k] = m;
  }
  slink->slink_seq2 = 0;
}

template<int verbose>
void ClusterSLINK<verbose>::update() {
  if (slink->Lambda[slink->slink_seq2] >= slink->M[slink->slink_seq2])
  {
    slink->M[slink->Pi[slink->slink_seq2]] = std::min(
        slink->M[slink->Pi[slink->slink_seq2]],
        slink->Lambda[slink->slink_seq2]);
    slink->Lambda[slink->slink_seq2] = slink->M[slink->slink_seq2];
    slink->Pi[slink->slink_seq2] = slink->slink_seq1;
  }
  else
  {
    slink->M[slink->Pi[slink->slink_seq2]] = std::min(
        slink->M[slink->Pi[slink->slink_seq2]],
        slink->M[slink->slink_seq2]);
  }
}

template <int verbose>
void ClusterSLINK<verbose>::finish_iter()
{
  while (slink->slink_seq2 < slink->slink_seq1)
  {
    update();
    slink->slink_seq2++;
  }
  for (j_t k = 0; k < slink->slink_seq1; k++)
  {
    if (slink->Lambda[k] >= slink->Lambda[slink->Pi[k]])
    {
      slink->Pi[k] = slink->slink_seq1;
    }
  }
}

template <int verbose>
void ClusterSLINK<verbose>::apply_ordered(j_t seq1, j_t seq2, d_t i)
{
  if (seq2 < 0 || seq2 >= n)
  {
    OPTIMOTU_CERR << "Sequence index" << seq2 << " out of range." << std::endl;
    OPTIMOTU_STOP("ClusterSLINK input error.");
  }
  if (seq1 < 0 || seq1 >= n)
  {
    OPTIMOTU_CERR << "Sequence index" << seq1 << " out of range." << std::endl;
    OPTIMOTU_STOP("ClusterSLINK input error.");
  }
  if (seq1 <= seq2)
  {
    OPTIMOTU_CERR << "ClusterSLINK requires seq2 < seq1.  seq2=" << seq2
                  << ", seq1=" << seq1 << std::endl;
    OPTIMOTU_STOP("ClusterSLINK input error.");
  }
  if (seq1 < slink->slink_seq1)
  {
    OPTIMOTU_CERR << "ClusterSLINK requires sequences in order."
                  << " Current seq1=" << seq1
                  << ", slink_seq1=" << slink->slink_seq1 << std::endl;
    OPTIMOTU_STOP("ClusterSLINK input error.");
  }
  while (seq1 > slink->slink_seq1)
  {
    finish_iter();
    ++slink->slink_seq1;
    slink->Pi[slink->slink_seq1] = slink->slink_seq1;
    init_iter();
  }
  if (seq2 < slink->slink_seq2)
  {
    OPTIMOTU_CERR << "ClusterSLINK requires sequences in order."
                  << " Current seq2=" << seq2
                  << ", slink_seq2=" << slink->slink_seq2 << std::endl;
    OPTIMOTU_STOP("ClusterSLINK input error.");
  }
  while (slink->slink_seq2 <= seq2)
  {
    if (slink->slink_seq2 == seq2)
    {
      slink->M[slink->slink_seq2] = std::min(i, slink->M[slink->slink_seq2]);
    }
    update();
    slink->slink_seq2++;
  }
}

template<int verbose>
void ClusterSLINK<verbose>::finalize() {
  if (uses_cache)
  {
    std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
    std::sort(
        merge_cache.begin(),
        merge_cache.end(),
        [](const MergeEdge &a, const MergeEdge &b)
        {
          if (a.seq1 != b.seq1)
          {
            return a.seq1 < b.seq1;
          }
          return a.seq2 < b.seq2;
        });
    // Tile MST pointer representations can emit the same parent pair more
    // than once; keep the lightest edge so replay stays in order.
    if (!merge_cache.empty())
    {
      std::size_t w = 0;
      for (std::size_t r = 1; r < merge_cache.size(); ++r)
      {
        if (merge_cache[r].seq1 == merge_cache[w].seq1 &&
            merge_cache[r].seq2 == merge_cache[w].seq2)
        {
          merge_cache[w].i = std::min(merge_cache[w].i, merge_cache[r].i);
        }
        else
        {
          merge_cache[++w] = merge_cache[r];
        }
      }
      merge_cache.resize(w + 1);
    }
    ensure_slink();
    for (const MergeEdge &e : merge_cache)
    {
      apply_ordered(e.seq1, e.seq2, e.i);
    }
    merge_cache.clear();
    merge_cache.shrink_to_fit();
    uses_cache = false;
    finish_iter();
    return;
  }
  ensure_slink();
  finish_iter();
}

template<int verbose>
ClusterSLINK<verbose> * ClusterSLINK<verbose>::make_inner_child(ClusterAlgorithm * parent, const j_t n) {
  // this is not locked, because it is called during the MappedClusterAlgorithm,
  // constructor, which is called by make_child(), which is already locked.
  this->uses_cache = true;
  auto child_ptr = new ClusterSLINK<verbose>(parent, n);
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  this->children.push_back(std::move(child));
  return child_ptr;
}

template<int verbose>
ClusterSLINK<verbose> * ClusterSLINK<verbose>::make_child() {
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  this->uses_cache = true;
  auto child_ptr = new ClusterSLINK<verbose>(this, n);
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  this->children.push_back(std::move(child));
  return child_ptr;
}

template<int verbose>
MappedClusterAlgorithm * ClusterSLINK<verbose>::make_child(PairGenerator * pg) {
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  this->uses_cache = true;
  auto child_ptr = new MappedClusterAlgorithmImpl<verbose>(this, pg);
  auto child = std::unique_ptr<ClusterAlgorithm>(
      (ClusterAlgorithm *)child_ptr);
  this->children.push_back(std::move(child));
  return child_ptr;
}

template <int verbose>
ClusterAlgorithm *ClusterSLINK<verbose>::make_child(
    const std::vector<std::size_t> &fwd_map)
{
  if (fwd_map.size() < 2)
  {
    return nullptr;
  }
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  this->uses_cache = true;
  auto rev_map = invert_map(fwd_map);
  auto child_ptr = new MappedClusterAlgorithmImpl<verbose>(
      this, fwd_map, rev_map);
  auto child = std::unique_ptr<ClusterAlgorithm>(
      (ClusterAlgorithm *)child_ptr);
  this->children.push_back(std::move(child));
  return child_ptr;
}

template <int verbose>
void ClusterSLINK<verbose>::operator()(j_t seq1, j_t seq2, d_t i, int thread)
{
  if (uses_cache)
  {
    if (seq1 == seq2)
    {
      return;
    }
    if (seq1 < seq2)
    {
      std::swap(seq1, seq2);
    }
    std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
    merge_cache.push_back({seq1, seq2, i});
    return;
  }
  ensure_slink();
  apply_ordered(seq1, seq2, i);
}

template<int verbose>
void ClusterSLINK<verbose>::operator()(PairGenerator & pg, d_t i, int thread) {
  if (pg) this->operator()(pg.i0(), pg.j0(), i, thread);
}

template <int verbose>
void ClusterSLINK<verbose>::write_to_matrix(internal_matrix_t &out)
{
  ensure_slink();
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  j_t j;
  std::size_t k = 0;
  for (std::uint32_t i = 0; i < this->n; i++) {
    j = i;
    d_t i2 = 0;
    while (i2 < this->m) {
      d_t max = slink->Lambda[j];
      if (this->m < max) max = this->m;
      while (i2 < max) {
        out[k++] = j;
        i2++;
      }
      if (i2 < this->m) {
        j = slink->Pi[j];
      }
    }
  }
}

template <int verbose>
void ClusterSLINK<verbose>::prepare_output()
{
}

template <int verbose>
void ClusterSLINK<verbose>::write_threshold_row(d_t t, int *dest) const
{
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  for (j_t i = 0; i < this->n; ++i)
  {
    j_t j = i;
    // Lambda is initialized to m, so this stops at the root of the
    // pointer representation. Same cut as write_to_matrix().
    while (slink->Lambda[j] <= t)
    {
      j = slink->Pi[j];
    }
    dest[i] = static_cast<int>(j);
  }
}

#ifdef OPTIMOTU_R

struct RevOrderElement {
  RevOrderElement * prev;
  int i;
};

template <int verbose>
Rcpp::List ClusterSLINK<verbose>::as_hclust(const Rcpp::CharacterVector &seqnames) const
{
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  Rcpp::IntegerMatrix merge(this->n - 1, 2);
  Rcpp::NumericVector height(this->n-1);
  Rcpp::IntegerVector order(this->n);

  // keep track of which cluster the active tips are in at each stage of clustering
  std::vector<int> clust_id;
  // ordering of earlier tips in the same cluster
  std::vector<RevOrderElement> ordering;
  // first element of each cluster (reverse_list does not know this)
  std::vector<RevOrderElement*> first;
  clust_id.reserve(this->n);
  ordering.reserve(this->n);
  first.reserve(this->n);

  std::vector<std::tuple<d_t, j_t, j_t>> slink_merge;
  slink_merge.reserve(this->n);
  int last_clust = 0;
  {
    for (int i = 0; i < (int)n;) {
      slink_merge.emplace_back(slink->Lambda[i], -i, -slink->Pi[i]);
      ++i;
      clust_id.push_back(-i);
      ordering.push_back({NULL, i});
      first.push_back(&ordering.back());
    }
    std::sort(slink_merge.begin(), slink_merge.end());

    for (const auto &c : slink_merge)
    {
      double d = this->dconv.inverse(std::get<0>(c));
      j_t i = -std::get<1>(c);
      j_t j = -std::get<2>(c);
      merge(last_clust, 0) = clust_id[i];
      merge(last_clust, 1) = clust_id[j];
      height[last_clust] = d;
      first[j]->prev = &ordering[i];
      first[j] = first[i];
      clust_id[j] = ++last_clust;
    }
  }
  int j = 0;
  for (j_t i = 0; i < this->n - 1; ++i) {
    if (i + (int)clust_id[i] == 1) {
      merge(last_clust, 0) = clust_id[i];
      merge(last_clust, 1) = clust_id[this->n - 1];
      height[last_clust] = 1.0;
      clust_id[this->n - 1] = ++last_clust;

      RevOrderElement * e = &ordering[i];
      while (e != NULL) {
        order[j] = e->i;
        ++j;
        e = e->prev;
      }
    }
  }

  Rcpp::List out;
  out["merge"] = merge;
  out["height"] = height;
  out["order"] = order;
  out["labels"] = seqnames;
  out["method"] = "single";
  out.attr("class") = "hclust";
  return out;
}
#endif // OPTIMOTU_R

template <int verbose>
void ClusterSLINK<verbose>::merge_into(DistanceConsumer &consumer)
{
  ensure_slink();
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  for (j_t j = 0; j < this->n; j++) {
    if (slink->Pi[j] != j)
    {
      consumer(slink->Pi[j], j, dconv.inverse(slink->Lambda[j]));
    }
  }
}

template<int verbose>
void ClusterSLINK<verbose>::merge_into(ClusterAlgorithm &consumer) {
  if (!consumer.accepts_unordered_pairs()) {
    OPTIMOTU_STOP(
      "ClusterSLINK::merge_into requires a consumer that accepts unordered pairs"
    );
  }
  ensure_slink();
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  for (j_t j = 0; j < this->n; j++) {
    if (slink->Pi[j] != j)
      consumer(slink->Pi[j], j, slink->Lambda[j]);
  }
}

template<int verbose>
void ClusterSLINK<verbose>::merge_into_parent() {
  if (!parent)
    return;
  this->merge_into(*parent);
}

template<int verbose>
double ClusterSLINK<verbose>::max_relevant(j_t seq1, j_t seq2, int thread) const {
  if (uses_cache && !slink)
  {
    return dconv.inverse(this->m - 1);
  }
  if (seq1 <= seq2) {
    OPTIMOTU_CERR << "ClusterSLINK requires seq2 < seq1.  seq2=" << seq2
                  << ", seq1=" << seq1 << std::endl;
    OPTIMOTU_STOP("ClusterSLINK input error.");
  }
  if (!slink)
  {
    return dconv.inverse(this->m - 1);
  }
  if (seq1 < slink->slink_seq1)
  {
    OPTIMOTU_CERR << "ClusterSLINK requires sequences in order."
                  << " Current seq1=" << seq1
                  << ", slink_seq1=" << slink->slink_seq1 << std::endl;
    OPTIMOTU_STOP("ClusterSLINK input error.");
  }
  if (seq1 > slink->slink_seq1)
  {
    return dconv.inverse(this->m - 1);
  }
  return dconv.inverse(slink->M[seq2] - 1);
}

template<int verbose>
double ClusterSLINK<verbose>::max_relevant(PairGenerator & pg, int thread) const {
  if (!pg) return dconv.inverse(this->m - 1);
  return this->max_relevant(pg.i0(), pg.j0(), thread);
}

template<int verbose>
std::unique_ptr<SingleClusterAlgorithm> create_cluster_slink(
  const DistanceConverter & dconv, j_t n)
{
  return std::make_unique<ClusterSLINK<verbose>>(dconv, n);
}

template<int verbose>
std::unique_ptr<SingleClusterAlgorithm> create_cluster_slink(
  const DistanceConverter & dconv, init_matrix_t & im)
{
  return std::make_unique<ClusterSLINK<verbose>>(dconv, im);
}

std::unique_ptr<SingleClusterAlgorithm> create_cluster_slink(
  const DistanceConverter & dconv, j_t n, int verbose)
{
  int v = (verbose > 4) ? 4 : verbose;
  if (v == 0) return create_cluster_slink<0>(dconv, n);
  if (v == 1) return create_cluster_slink<1>(dconv, n);
  if (v == 2) return create_cluster_slink<2>(dconv, n);
  if (v == 3) return create_cluster_slink<3>(dconv, n);
  return create_cluster_slink<4>(dconv, n);
}

std::unique_ptr<SingleClusterAlgorithm> create_cluster_slink(
  const DistanceConverter & dconv, init_matrix_t & im, int verbose)
{
  int v = (verbose > 4) ? 4 : verbose;
  if (v == 0) return create_cluster_slink<0>(dconv, im);
  if (v == 1) return create_cluster_slink<1>(dconv, im);
  if (v == 2) return create_cluster_slink<2>(dconv, im);
  if (v == 3) return create_cluster_slink<3>(dconv, im);
  return create_cluster_slink<4>(dconv, im);
}
