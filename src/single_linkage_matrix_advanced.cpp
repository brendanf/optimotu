#include <fstream>
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "ClusterMatrix.h"

Rcpp::IntegerMatrix single_linkage_matrix2(
    const std::string file,
    const Rcpp::CharacterVector &seqnames,
    const DistanceConverter &dconv,
    const int m,
    const int do_binary_search = 3
) {
  Rcpp::IntegerMatrix im(m, seqnames.size());
  ClusterMatrix<RcppParallel::RMatrix<int>, Rcpp::IntegerMatrix> cm(dconv, im, do_binary_search);
  std::ifstream infile(file);
  std::size_t seq1, seq2;
  double dist;
  while(infile >> seq1 >> seq2 >> dist) {
    cm(seq1, seq2, dist);
  }
  return im;
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix single_linkage_matrix2_uniform(
      const std::string file,
      const Rcpp::CharacterVector &seqnames,
      const float dmin,
      const float dmax,
      const float dstep,
      const uint8_t do_binary_search = 3
) {
   const UniformDistanceConverter dconv(dmin, dstep);
   const int m = (int) ceilf((dmax - dmin)/dstep) + 1;
   return single_linkage_matrix2(file, seqnames, dconv, m, do_binary_search);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix single_linkage_matrix2_array(
      const std::string file,
      const Rcpp::CharacterVector &seqnames,
      const std::vector<double> &thresholds,
      const uint8_t do_binary_search = 3
) {
   const ArrayDistanceConverter dconv(thresholds);
   const int m = thresholds.size();
   return single_linkage_matrix2(file, seqnames, dconv, m, do_binary_search);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix single_linkage_matrix2_cached(
      const std::string file,
      const Rcpp::CharacterVector &seqnames,
      const std::vector<double> &thresholds,
      const double precision,
      const uint8_t do_binary_search = 3
) {
   const CachedDistanceConverter dconv(thresholds, precision);
   const int m = thresholds.size();
   return single_linkage_matrix2(file, seqnames, dconv, m, do_binary_search);
}
