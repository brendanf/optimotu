// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include <deque>
#include <array>
#include <cstring>
#include "ClusterMatrix.h"
#include "optimotu.h"

template<bool BM, int F, typename A, int V>
ClusterMatrix<BM, F, A, V>::ClusterMatrix(ClusterAlgorithm * parent, const j_t n):
  SingleClusterAlgorithm(parent, n),
  clust_array(m*n),
  ca(&clust_array[0]),
  toclust(m, 0)
{
  initialize();
}

template<bool BM, int F, typename A, int V>
ClusterMatrix<BM, F, A, V>::ClusterMatrix(const DistanceConverter &dconv, size_t n) :
  SingleClusterAlgorithm(dconv, n),
  clust_array(m*n), ca(&clust_array[0]), toclust(m, 0) {
  initialize();
}

template<bool BM, int F, typename A, int V>
ClusterMatrix<BM, F, A, V>::ClusterMatrix(
    const DistanceConverter &dconv,
    init_matrix_t im
) :
  SingleClusterAlgorithm(dconv, im),
  clust_array(im), ca(&clust_array[0]), toclust(m, 0) {
  initialize();
}

template<bool BM, int F, typename A, int V>
void ClusterMatrix<BM, F, A, V>::initialize() {
  auto ca = clust_array.begin();
  for (j_t j = 0; j < n; j++) {
    for (d_t i = 0; i < m; i++) {
      *ca = j;
      ++ca;
    }
  }
}

template<bool BM, int F, typename A, int V>
void ClusterMatrix<BM, F, A, V>::operator()(j_t seq1, j_t seq2, d_t i, int thread) {
  if (i >= m) return;
  if (seq1 == seq2) return;

  if (seq2 < 0 || seq2 >= n) {
    OPTIMOTU_CERR << "Sequence index" << seq2 << " out of range." << std::endl;
    OPTIMOTU_STOP("ClusterMatrix input error.");
  }
  if (seq1 < 0 || seq1 >= n) {
    OPTIMOTU_CERR << "Sequence index" << seq1 << " out of range." << std::endl;
    OPTIMOTU_STOP("ClusterMatrix input error.");
  }

  {
    std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
    if (clust_array[i + seq1*m] == clust_array[i + seq2*m]) return;
  }
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  // check again in case the lock had to be released
  if (clust_array[i + seq1*m] == clust_array[i + seq2*m]) return;

  // imin == last row index where seq1 and seq2 are in separate clusters
  // imax == first row index where seq1 and seq2 are in the same cluster
  d_t imin = i, imax = m;
  // shortcut for case when the clusters are not yet joined at any level
  if (clust_array[m-1 + seq1*m] != clust_array[m-1 + seq2*m]) {
    imin = m-1;
  }
  int *c1min = ca + imin + seq1*m;
  int *c2min = ca + imin + seq2*m;
  int *c1max = ca + imax + seq1*m;
  if (BM) {
    // binary search to find imax and imin
    while (c1max - c1min > 1) {
      int diff = (c1max-c1min) / 2;
      int *c1mean = c1min + diff;
      int *c2mean = c2min + diff;
      if (*c1mean == *c2mean) {
        c1max=c1mean;
      } else {
        c1min=c1mean;
        c2min=c2mean;
      }
    }
    imax=c1max - (ca + seq1*m);
    imin=imax-1;
  } else {
    // OPTIMOTU_COUT << "initializing linear search for join of " << seq1 << " and "
    //             << seq2 << " at i=" << i << std::endl;
    auto c1 = ca + std::min(seq1, seq2)*m + i;
    auto c2 = ca + std::max(seq1, seq2)*m + i;
    // OPTIMOTU_COUT << "searching";
    imax=i;
    while(imax < m && *c2 != *c1) {
      // OPTIMOTU_COUT << ".";
      ++c1;
      ++c2;
      ++imax;
    }
    // OPTIMOTU_COUT << "done" << std::endl;
    imin = imax - 1;
  }
  // OPTIMOTU_COUT << " imax=" << imax << std::endl;
  size_t j1m, j2m;
  if (seq1 > seq2) {
    j1m = seq2*m;
    j2m = seq1*m;
  } else {
    j1m = seq1*m;
    j2m = seq2*m;
  }
  // cache the new clusters for seq1 and seq2 at each level where they
  // need to be joined
  // OPTIMOTU_COUT << "newclust=";
  for (d_t ii = i; ii < imax; ii++) {
    toclust[ii] = std::min(clust_array[ii + j1m], clust_array[ii + j2m]);
    // OPTIMOTU_COUT << toclust[ii] << " ";
  }
  // OPTIMOTU_COUT << "done" << std::endl;
  // for each sequence which is clustered with seq1 or seq2 at imin,
  // assign it to newclust in the same range where it is
  std::size_t jjmin = std::max(clust_array[j1m + imin], clust_array[j2m + imin])*m;
  std::size_t jjmax = n*m;
  for (std::size_t jj = jjmin; jj < jjmax; jj += m) {
    if (jj == j1m) continue;
    if (jj == j2m) continue;
    // save j1m and j2m for last because we use them for comparisons to others
    size_t jcomp;
    auto ii = ca + jj + imin;
    if (*ii == clust_array[j1m + imin]) {
      jcomp = j1m;
    } else if (*ii == clust_array[j2m + imin]) {
      jcomp = j2m;
    } else {
      continue;
    }
    if (F==TOPDOWN_FILL) {
      int *icomp = ca + jcomp + imin;
      int *ito = &toclust[imin];
      while (*ii == *icomp) {
        *ii = *ito;
        --ii;
        --ito;
        --icomp;
        if (ito < &toclust[i]) break;
      }
    } else {
      d_t copystart;
      auto ca = clust_array.begin() + jj;
      // OPTIMOTU_COUT << " setting cluster " << (jj / m);
      if (F==BINARY_FILL) {
        // iimin == last value of ii where this seq is not clustered with seq1 or seq2
        // copystart == first value of ii where this seq is clustered with seq1 or seq2
        size_t iimin = i;
        copystart = imin;
        // if it is already clustered at i, then shortcut.
        if (clust_array[jj + i] == clust_array[jcomp + i]) {
          copystart = i;
        }
        // binary search to find iimin and iimax
        while (copystart - iimin > 1) {
          size_t iimean = (copystart + iimin)/2;
          if (clust_array[jj + iimean] == clust_array[jcomp + iimean]) {
            copystart = iimean;
          } else {
            iimin = iimean;
          }
        }
        ca += copystart;
        // OPTIMOTU_COUT << " from " << iimax << " to " << imin << std::endl;
        // assign
      } else {
        // OPTIMOTU_COUT << "starting fill linear search for j=" << jj/m << " i=" << i << std::endl;
        copystart = i;
        ca += i;
        auto jc = clust_array.begin() + jcomp + i;
        while (*ca != *jc) {
          // OPTIMOTU_COUT << "ca=" << *ca << " jc=" << *jc << " j=" << copystart << std::endl;
          ++ca;
          ++jc;
          ++copystart;
        }
        // OPTIMOTU_COUT << "cluster " << jj/m << " matching to " << jcomp/m <<
        // "; updating from " << copystart << " to " << imin << std::endl;
      }
      std::memcpy(
        &*ca,
        &*(toclust.begin()) + copystart,
        (imax - copystart) * sizeof(int)
      );
    }
  }
  std::memcpy(
    &*clust_array.begin() + j1m + i,
    &*toclust.begin() + i,
    (imax - i) * sizeof(int)
  );
  std::memcpy(
    &*clust_array.begin() + j2m + i,
    &*toclust.begin() + i,
    (imax - i) * sizeof(int)
  );
}

template<bool BM, int F, typename A, int V>
void ClusterMatrix<BM, F, A, V>::merge_into(DistanceConsumer &consumer) {
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  for (size_t j = 1; j < n; ++j) {
    j_t i = j*m;
    if (clust_array[i + m - 1] == (int)j) continue;
    while(clust_array[i] == (int)j) {
      ++i;
    }
    consumer(j, clust_array[i], dconv.inverse(i - j*m));
  }
}

template<bool BM, int F, typename A, int V>
void ClusterMatrix<BM, F, A, V>::merge_into(ClusterAlgorithm &consumer) {
  if (!consumer.accepts_unordered_pairs()) {
    OPTIMOTU_STOP(
      "ClusterMatrix::merge_into requires a consumer that accepts unordered pairs"
    );
  }
  // TODO check that the distance converters are really compatible
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  for (j_t j = 1; j < n; ++j) {
    j_t i = j*m;
    if (clust_array[i + m - 1] == (int)j) continue;
    while(clust_array[i] == (int)j) {
      ++i;
    }
    consumer(j, clust_array[i], (d_t)(i - j*m));
  }
}

template<bool BM, int F, typename A, int V>
SingleClusterAlgorithm * ClusterMatrix<BM, F, A, V>::make_child() {
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto child_ptr = new ClusterMatrix<BM, F, A, V>(this, n);
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  this->children.push_back(std::move(child));
  return child_ptr;
}

#define CIM_REF_SPEC(BM, F) \
template<> SingleClusterAlgorithm * ClusterMatrix<BM, F, internal_matrix_ref_t, 0>::make_inner_child(ClusterAlgorithm * parent, const j_t n) { \
  /* Called while make_child() already holds this mutex. */ \
  auto child_ptr = new ClusterMatrix<BM, F, std::vector<int>, 0>(parent, n); \
  auto child = std::unique_ptr<ClusterAlgorithm>( (ClusterAlgorithm*)child_ptr ); \
  this->children.push_back(std::move(child)); \
  return child_ptr; \
} \
template<> SingleClusterAlgorithm * ClusterMatrix<BM, F, internal_matrix_ref_t, 0>::make_child() { \
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex); \
  auto child_ptr = new ClusterMatrix<BM, F, std::vector<int>, 0>(this, n); \
  auto child = std::unique_ptr<ClusterAlgorithm>( (ClusterAlgorithm*)child_ptr ); \
  this->children.push_back(std::move(child)); \
  return child_ptr; \
}
CIM_REF_SPEC(true, LINEAR_FILL)
CIM_REF_SPEC(true, BINARY_FILL)
CIM_REF_SPEC(true, TOPDOWN_FILL)
CIM_REF_SPEC(false, LINEAR_FILL)
CIM_REF_SPEC(false, BINARY_FILL)
CIM_REF_SPEC(false, TOPDOWN_FILL)
#undef CIM_REF_SPEC

template<bool BM, int F, typename A, int V>
MappedClusterAlgorithm * ClusterMatrix<BM, F, A, V>::make_child(
  PairGenerator * pg
) {
  // std::lock_guard<std::mutex> lock(this->mutex);
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto child_ptr = new MappedClusterAlgorithmImpl<V>(this, pg);
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  this->children.push_back(std::move(child));
  return child_ptr;
}

template <bool BM, int F, typename A, int V>
ClusterAlgorithm *ClusterMatrix<BM, F, A, V>::make_child(
    const std::vector<std::size_t> &fwd_map)
{
  if (fwd_map.size() < 2)
  {
    return nullptr;
  }
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  auto rev_map = invert_map(fwd_map);
  auto child_ptr = new MappedClusterAlgorithmImpl<V>(this, fwd_map, rev_map);
  auto child = std::unique_ptr<ClusterAlgorithm>(
      (ClusterAlgorithm *)child_ptr);
  this->children.push_back(std::move(child));
  return child_ptr;
}

template<bool BM, int F, typename A, int V>
SingleClusterAlgorithm * ClusterMatrix<BM, F, A, V>::make_inner_child(ClusterAlgorithm * parent, const j_t n) {
  // This is called from MappedClusterAlgorithm construction while make_child()
  // already holds this object's mutex. Locking again here can deadlock.
  auto child_ptr = new ClusterMatrix<BM, F, A, V>(parent, n);
  auto child = std::unique_ptr<ClusterAlgorithm>(
    (ClusterAlgorithm*)child_ptr
  );
  this->children.push_back(std::move(child));
  return child_ptr;
}

template<bool BM, int F, typename A, int V>
double ClusterMatrix<BM, F, A, V>::max_relevant(j_t seq1, j_t seq2, int thread) const {
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  if (seq1 == seq2) return 0.0;
  j_t j1 = seq1*m, j2 = seq2*m;
  auto c1 = ca + j1, c2 = ca + j2,
    c1max = c1+m;
  if (clust_array[j1 + m - 1] != clust_array[j2 + m - 1]) {
    return dconv.inverse(m-1);
  }
  if (*c1 == *c2) return -1;
  if (BM) {
    // binary search to find imax and imin
    while (c1max - c1 > 1) {
      int diff = (c1max-c1) / 2;
      auto c1mean = c1 + diff;
      auto c2mean = c2 + diff;
      if (*c1mean == *c2mean) {
        c1max=c1mean;
      } else {
        c1=c1mean;
        c2=c2mean;
      }
    }
    return dconv.inverse(c1 - ca - j1);
  } else {
    // linear search version
    while (c1 < c1max && *c1 != *c2) {
      ++c1;
      ++c2;
    }
    return dconv.inverse(c1 - ca - j1 - 1);
  }
}

template<bool BM, int F, typename A, int V>
void ClusterMatrix<BM, F, A, V>::write_to_matrix(internal_matrix_t &out) {
  // shortcut if this is already our data matrix
  if (intptr_t(&out[0]) == intptr_t(&clust_array[0])) return;
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  std::copy(clust_array.begin(), clust_array.end(), out.begin());
}

#ifdef OPTIMOTU_R

struct FwdOrderElement {
  FwdOrderElement * next;
  int i;
};

template<bool BM, int F, typename A, int V>
Rcpp::List ClusterMatrix<BM, F, A, V>::as_hclust(
    const Rcpp::CharacterVector &seqnames
) const {
  // output objects
  Rcpp::IntegerMatrix merge(this->n - 1, 2);
  Rcpp::NumericVector height(this->n-1);
  Rcpp::IntegerVector order(this->n);

  // keep track of which cluster the active tips are in at each stage of clustering
  std::vector<int> clust_id;
  // ordering of later tips in the same cluster
  std::vector<FwdOrderElement> ordering;
  // last element of each cluster (forward_list does not know this)
  std::vector<FwdOrderElement*> last;
  clust_id.reserve(this->n);
  ordering.reserve(this->n);
  last.reserve(this->n);
  // which tips are still "active"
  std::array<std::deque<int>, 2> remaining;
  // keep two versions
  int this_remaining = 0, next_remaining = 1;
  for (int i = 0; i < (int)n;) {
    remaining[this_remaining].push_back(i);
    ++i;
    clust_id.push_back(-i);
    ordering.push_back({NULL, i});
    last.push_back(&ordering.back());
  }
  int last_clust = 0;
  {
    std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
    for (j_t j = 0; j < this->m * this->n; j += this->n) {
      double d = this->dconv.inverse(j);
      for (int i : remaining[this_remaining]) {
        int clust = this->ca[j+i];
        if (clust == i) {
          remaining[next_remaining].push_back(i);
        } else {
          merge(last_clust, 0) = clust_id[i];
          merge(last_clust, 1) = clust_id[clust];
          height[last_clust] = d;
          last[clust]->next = &ordering[i];
          last[clust] = last[i];
          clust_id[clust] = ++last_clust;
        }
        remaining[this_remaining].pop_front();
      }
      next_remaining++;
      next_remaining %= 2;
      this_remaining++;
      this_remaining %= 2;
    }
  }
  int j = 0;
  for (int i : remaining[this_remaining]) {
    if (i > 0) {
      merge(last_clust, 0) = i;
      merge(last_clust, 1) = clust_id[0];
      height[last_clust] = 1.0;
      clust_id[0] = ++last_clust;
    }
    FwdOrderElement * e = &ordering[i];
    while (e != NULL) {
      order[j] = e->i;
      ++j;
      e = e->next;
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



template<>
ClusterMatrix<true, LINEAR_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<true, BINARY_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<true, TOPDOWN_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<false, LINEAR_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<false, BINARY_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<false, TOPDOWN_FILL>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<true, LINEAR_FILL, std::vector<int>, 1>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<true, BINARY_FILL, std::vector<int>, 1>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<true, TOPDOWN_FILL, std::vector<int>, 1>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<false, LINEAR_FILL, std::vector<int>, 1>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<false, BINARY_FILL, std::vector<int>, 1>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<false, TOPDOWN_FILL, std::vector<int>, 1>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<true, LINEAR_FILL, std::vector<int>, 2>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<true, BINARY_FILL, std::vector<int>, 2>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<true, TOPDOWN_FILL, std::vector<int>, 2>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<false, LINEAR_FILL, std::vector<int>, 2>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<false, BINARY_FILL, std::vector<int>, 2>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<false, TOPDOWN_FILL, std::vector<int>, 2>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<true, LINEAR_FILL, std::vector<int>, 3>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<true, BINARY_FILL, std::vector<int>, 3>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<true, TOPDOWN_FILL, std::vector<int>, 3>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<false, LINEAR_FILL, std::vector<int>, 3>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<false, BINARY_FILL, std::vector<int>, 3>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<false, TOPDOWN_FILL, std::vector<int>, 3>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

template<>
ClusterMatrix<true, LINEAR_FILL, std::vector<int>, 4>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<true, BINARY_FILL, std::vector<int>, 4>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<true, TOPDOWN_FILL, std::vector<int>, 4>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<false, LINEAR_FILL, std::vector<int>, 4>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<false, BINARY_FILL, std::vector<int>, 4>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;
template<>
ClusterMatrix<false, TOPDOWN_FILL, std::vector<int>, 4>::ClusterMatrix(
    const DistanceConverter &dconv, init_matrix_t
) = delete;

#define deleted_funcs template<>                                     \
cm::ClusterMatrix(ClusterAlgorithm * parent, const j_t n) = delete;         \
template<>                                                           \
cm::ClusterMatrix(const DistanceConverter &dconv, size_t n) = delete

#define cm ClusterMatrix<true, LINEAR_FILL, internal_matrix_ref_t, 0>
deleted_funcs;
#undef cm
#define cm ClusterMatrix<true, BINARY_FILL, internal_matrix_ref_t, 0>
deleted_funcs;
#undef cm
#define cm ClusterMatrix<true, TOPDOWN_FILL, internal_matrix_ref_t, 0>
deleted_funcs;
#undef cm
#define cm ClusterMatrix<false, LINEAR_FILL, internal_matrix_ref_t, 0>
deleted_funcs;
#undef cm
#define cm ClusterMatrix<false, BINARY_FILL, internal_matrix_ref_t, 0>
deleted_funcs;
#undef cm
#define cm ClusterMatrix<false, TOPDOWN_FILL, internal_matrix_ref_t, 0>
deleted_funcs;
#undef cm
#undef deleted_funcs

template<int verbose>
std::unique_ptr<SingleClusterAlgorithm> create_cluster_matrix(
  const DistanceConverter & dconv, j_t n, bool binary_search, int fill_type)
{
  if (binary_search) {
    switch (fill_type) {
    case LINEAR_FILL:
      return std::make_unique<ClusterMatrix<true, LINEAR_FILL, std::vector<int>, verbose>>(dconv, n);
    case BINARY_FILL:
      return std::make_unique<ClusterMatrix<true, BINARY_FILL, std::vector<int>, verbose>>(dconv, n);
    case TOPDOWN_FILL:
      return std::make_unique<ClusterMatrix<true, TOPDOWN_FILL, std::vector<int>, verbose>>(dconv, n);
    default:
      OPTIMOTU_STOP("unknown fill type");
    }
  } else {
    switch (fill_type) {
    case LINEAR_FILL:
      return std::make_unique<ClusterMatrix<false, LINEAR_FILL, std::vector<int>, verbose>>(dconv, n);
    case BINARY_FILL:
      return std::make_unique<ClusterMatrix<false, BINARY_FILL, std::vector<int>, verbose>>(dconv, n);
    case TOPDOWN_FILL:
      return std::make_unique<ClusterMatrix<false, TOPDOWN_FILL, std::vector<int>, verbose>>(dconv, n);
    default:
      OPTIMOTU_STOP("unknown fill type");
    }
  }
}

std::unique_ptr<SingleClusterAlgorithm> create_cluster_matrix(
  const DistanceConverter & dconv, j_t n,
  bool binary_search, int fill_type, int verbose)
{
  int v = (verbose > 4) ? 4 : verbose;
  if (v == 0) return create_cluster_matrix<0>(dconv, n, binary_search, fill_type);
  if (v == 1) return create_cluster_matrix<1>(dconv, n, binary_search, fill_type);
  if (v == 2) return create_cluster_matrix<2>(dconv, n, binary_search, fill_type);
  if (v == 3) return create_cluster_matrix<3>(dconv, n, binary_search, fill_type);
  return create_cluster_matrix<4>(dconv, n, binary_search, fill_type);
}

std::unique_ptr<SingleClusterAlgorithm> create_cluster_matrix(
  const DistanceConverter & dconv, init_matrix_t & im,
  bool binary_search, int fill_type)
{
  if (binary_search) {
    switch (fill_type) {
    case LINEAR_FILL:
      return std::make_unique<ClusterMatrix<true, LINEAR_FILL, internal_matrix_ref_t, 0>>(dconv, im);
    case BINARY_FILL:
      return std::make_unique<ClusterMatrix<true, BINARY_FILL, internal_matrix_ref_t, 0>>(dconv, im);
    case TOPDOWN_FILL:
      return std::make_unique<ClusterMatrix<true, TOPDOWN_FILL, internal_matrix_ref_t, 0>>(dconv, im);
    default:
      OPTIMOTU_STOP("unknown fill type");
    }
  } else {
    switch (fill_type) {
    case LINEAR_FILL:
      return std::make_unique<ClusterMatrix<false, LINEAR_FILL, internal_matrix_ref_t, 0>>(dconv, im);
    case BINARY_FILL:
      return std::make_unique<ClusterMatrix<false, BINARY_FILL, internal_matrix_ref_t, 0>>(dconv, im);
    case TOPDOWN_FILL:
      return std::make_unique<ClusterMatrix<false, TOPDOWN_FILL, internal_matrix_ref_t, 0>>(dconv, im);
    default:
      OPTIMOTU_STOP("unknown fill type");
    }
  }
}
