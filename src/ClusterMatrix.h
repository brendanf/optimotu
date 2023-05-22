#ifndef _ClusterMatrix_
#define _ClusterMatrix_

#include "ClusterAlgorithm.h"
#include <RcppParallel.h>

#define BINARY_MAX  (uint8_t)0b00000001
#define BINARY_FILL (uint8_t)0b00000010
#define TOPDOWN_FILL (uint8_t)0b00000100

// Matrix representation of a hierarchical clustering of n items at m different
// thresholds.
template <class ARRAY_T = std::vector<int>, class PARENT_T = std::vector<int>>
struct ClusterMatrix : public ClusterAlgorithm {
private:
  ARRAY_T clust_array;
  int *const ca;
  std::vector<int> toclust;
  // std::vector<int> fromclust;
  const uint8_t do_binary_search;

public:
  ClusterMatrix(const DistanceConverter &dconv, size_t n, int m,
                uint8_t do_binary_search = 3) :
  ClusterAlgorithm(dconv, n, m), clust_array(m*n), ca(&clust_array[0]), toclust(m, 0),
  do_binary_search(do_binary_search) {
    auto ca = clust_array.begin();
    for (j_t j = 0; j < n; j++) {
      for (d_t i = 0; i < m; i++) {
        *ca = j;
        ++ca;
      }
    }
  };
  ClusterMatrix(const DistanceConverter &dconv, PARENT_T &parent,
                uint8_t do_binary_search = 3);

  using ClusterAlgorithm::operator();

  void operator()(j_t seq1, j_t seq2, d_t i) override {
    if (i >= m) return;
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
    int *c2max = ca + imax + seq2*m;
    if (do_binary_search & BINARY_MAX) {
      // binary search to find imax and imin
      while (c1max - c1min > 1) {
        int diff = (c1max-c1min) / 2;
        int *c1mean = c1min + diff;
        int *c2mean = c2min + diff;
        if (*c1mean == *c2mean) {
          c1max=c1mean;
          c2max=c2mean;
        } else {
          c1min=c1mean;
          c2min=c2mean;
        }
      }
      imax=c1max - (ca + seq1*m);
      imin=imax-1;
    } else {
      // Rcpp::Rcout << "initializing linear search for join of " << seq1 << " and "
      //             << seq2 << " at i=" << i << std::endl;
      auto c1 = ca + std::min(seq1, seq2)*m + i;
      auto c2 = ca + std::max(seq1, seq2)*m + i;
      // Rcpp::Rcout << "searching";
      imax=i;
      while(imax < m && *c2 != *c1) {
        // Rcpp::Rcout << ".";
        ++c1;
        ++c2;
        ++imax;
      }
      // Rcpp::Rcout << "done" << std::endl;
      imin = imax - 1;
    }
    // Rcpp::Rcout << " imax=" << imax << std::endl;
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
    // Rcpp::Rcout << "newclust=";
    for (d_t ii = i; ii < imax; ii++) {
      toclust[ii] = std::min(clust_array[ii + j1m], clust_array[ii + j2m]);
      // Rcpp::Rcout << toclust[ii] << " ";
    }
    // Rcpp::Rcout << "done" << std::endl;
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
      if (do_binary_search & TOPDOWN_FILL) {
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
      // Rcpp::Rcout << " setting cluster " << (jj / m);
      if (do_binary_search & BINARY_FILL) {
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
      // Rcpp::Rcout << " from " << iimax << " to " << imin << std::endl;
      // assign
      } else {
        // Rcpp::Rcout << "starting fill linear search for j=" << jj/m << " i=" << i << std::endl;
        copystart = i;
        ca += i;
        auto jc = clust_array.begin() + jcomp + i;
        while (*ca != *jc) {
          // Rcpp::Rcout << "ca=" << *ca << " jc=" << *jc << " j=" << copystart << std::endl;
          ++ca;
          ++jc;
          ++copystart;
        }
        // Rcpp::Rcout << "cluster " << jj/m << " matching to " << jcomp/m <<
          // "; updating from " << copystart << " to " << imin << std::endl;
      }
      memcpy(
        &*ca,
        &*(toclust.begin()) + copystart,
        (imax - copystart) * sizeof(int)
      );
    }
    }
    memcpy(
      &*clust_array.begin() + j1m + i,
      &*toclust.begin() + i,
      (imax - i) * sizeof(int)
    );
    memcpy(
      &*clust_array.begin() + j2m + i,
      &*toclust.begin() + i,
      (imax - i) * sizeof(int)
    );
  };

  void merge_into(DistanceConsumer &consumer) const override {
    for (j_t j = 1; j < n; ++j) {
      j_t i = j*m, imax = j*(m+1);
      while(clust_array[i] == j && i < imax) {
        ++i;
      }
      if (i < imax) consumer(j, clust_array[i], dconv.inverse(i - j*m));
    }
  };

  void merge_into(ClusterAlgorithm &consumer) const override {
    // TODO check that the distance converters are really compatible
    for (j_t j = 1; j < n; ++j) {
      j_t i = j*m, imax = j*(m+1);
      while(clust_array[i] == j && i < imax) {
        ++i;
      }
      if (i < imax) consumer(j, clust_array[i], (d_t)(i - j*m));
    }
  };

  double max_relevant(j_t seq1, j_t seq2) const override {
    if (seq1 == seq2) return 0.0;
    j_t j1 = seq1*m, j2 = seq2*m;
    auto c1 = clust_array.begin() + j1, c2 = clust_array.begin() + j2,
      c1max = c1+m;
    if (clust_array[j1 + m - 1] == clust_array[j2 + m - 1]) {
      return dconv.inverse(m-1);
    }
    // linear search version
    while (c1 < c1max && *c1 == *c2) {
      ++c1;
      ++c2;
    }
    // maybe try binary search?
    return dconv.inverse(c1 - clust_array.begin() - j1 - 1);
  };
};

template<>
ClusterMatrix<RcppParallel::RMatrix<int>, Rcpp::IntegerMatrix>::ClusterMatrix(
  const DistanceConverter &dconv, Rcpp::IntegerMatrix &im,
  uint8_t do_binary_search
) :
  ClusterAlgorithm(dconv, im.ncol(), im.nrow()), clust_array(im), ca(&clust_array[0]),
  toclust(m, 0), do_binary_search(do_binary_search) {
  auto ca = clust_array.begin();
  for (j_t j = 0; j < n; j++) {
    for (d_t i = 0; i < m; i++) {
      *ca = j;
      ++ca;
    }
  }
};

#endif
