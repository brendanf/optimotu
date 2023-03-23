#include <fstream>
#include <Rcpp.h>
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "DistanceConverter.h"

struct MergeWorker
{
  RcppParallel::RMatrix<int> clust_array;
  const std::size_t m, n;
  std::size_t j1m, j2m;

  MergeWorker(Rcpp::IntegerMatrix clust_array,
              const std::size_t m,
              const std::size_t n,
              const std::size_t j1,
              const std::size_t j2)
    : clust_array(clust_array), m(m), n(n), j1m(j1*m), j2m(j2*m) {};

  void operator()(std::size_t i) {
    int from = clust_array[i + j2m];
    int to = clust_array[i + j1m];
    if (to > from) {
      from = to;
      to = clust_array[i + j2m];
    }
    int xmax = n*m;
    for (int x = i + from*m; x < xmax; x+=m) {
      if (clust_array[x] == from) clust_array[x] = to;
    }
  }
};


Rcpp::IntegerMatrix single_linkage_matrix(
    const std::string file,
    const Rcpp::CharacterVector &seqnames,
    const DistanceConverter &dconv,
    const int m,
    const int threads=1,
    const int minsplit=1
) {
  RcppThread::ThreadPool::globalInstance().setNumThreads(threads);
  const std::size_t n = seqnames.size();
  // pointer to which cluster a sequence belongs to at each threshold
  Rcpp::IntegerMatrix output(m, n);
  RcppParallel::RMatrix<int> clust_array(output);
  std::size_t i, j, j1m, j2m;
  for (i = 0; i < (size_t)m; i++) {
    for (j = 0; j < n; j++) {
      clust_array[i + j*m] = j;
    }
  }
  std::ifstream infile(file);
  std::size_t seq1, seq2;
  double dist;
  while(infile >> seq1 >> seq2 >> dist) {
    if (seq1 == seq2) continue;
    i = dconv.convert(dist);
    if (i >= (size_t)m) continue;
    if (clust_array[i + seq1*m] == clust_array[i + seq2*m]) continue;
    std::size_t imin = i, imax = m;
    // shortcut for case when the clusters are not yet joined at any level
    if (clust_array[m-1 + seq1*m] != clust_array[m-1 + seq2*m]) {
      imin = m-1;
    }
    while (imax - imin > 1) {
      std::size_t imean = (imin + imax) / 2;
      if (clust_array[imean + seq1*m] == clust_array[imean + seq2*m]) {
        imax = imean;
      } else {
        imin = imean;
      }
    }
    if (seq1 > seq2) {
      j1m = seq2*m;
      j2m = seq1*m;
    } else {
      j1m = seq1*m;
      j2m = seq2*m;
    }
    if (threads > 1) {
    RcppThread::parallelFor(i, imax, [&clust_array, j1m, j2m, n, m] (std::size_t i) {
      int from = clust_array[i + j2m];
      int to = clust_array[i + j1m];
      if (to > from) {
        from = to;
        to = clust_array[i + j2m];
      }
      int xmax = n*m;
      for (int x = i + from*m; x < xmax; x+=m) {
        if (clust_array[x] == from) clust_array[x] = to;
      }
    });
    } else {
      for (std::size_t ii = i; ii < imax; ii++) {
        int from = clust_array[ii + j2m];
        int to = clust_array[ii + j1m];
        if (to > from) {
          from = to;
          to = clust_array[ii + j2m];
        }
        int xmax = n*m;
        for (int x = ii + from*m; x < xmax; x+=m) {
          if (clust_array[x] == from) clust_array[x] = to;
        }
      }
    }
  }
  return output;
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix single_linkage_matrix_uniform(
      const std::string file,
      const Rcpp::CharacterVector &seqnames,
      const float dmin,
      const float dmax,
      const float dstep,
      const int threads=1,
      const int minsplit=1
) {
   const UniformDistanceConverter dconv(dmin, dstep);
   const int m = (int) ceilf((dmax - dmin)/dstep) + 1;
   return single_linkage_matrix(file, seqnames, dconv, m, threads, minsplit);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix single_linkage_matrix_array(
      const std::string file,
      const Rcpp::CharacterVector &seqnames,
      const std::vector<double> &thresholds,
      const int threads=1,
      const int minsplit=1
) {
   const ArrayDistanceConverter dconv(thresholds);
   const int m = thresholds.size();
   return single_linkage_matrix(file, seqnames, dconv, m, threads, minsplit);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix single_linkage_matrix_cached(
      const std::string file,
      const Rcpp::CharacterVector &seqnames,
      const std::vector<double> &thresholds,
      const double precision,
      const int threads=1,
      const int minsplit=1
) {
   const CachedDistanceConverter dconv(thresholds, precision);
   const int m = thresholds.size();
   return single_linkage_matrix(file, seqnames, dconv, m, threads, minsplit);
}
