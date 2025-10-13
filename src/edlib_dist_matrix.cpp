#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>

#include "pairwise_alignment.h"
#include "SparseDistanceMatrix.h"
#include "EdlibDistWorker.h"

#include <cstdint>

//' Sparse distance matrix between DNA sequences
//'
//' @name seq_distmx
//'
//' @param seq (`character` vector) DNA sequences to calculate distances for
//' @param dist_threshold (`numeric` scalar) maximum sequence distance (edit
//' distance / alignment length) threshold for reporting
//' @param constrain (`logical` flag) if `TRUE`, the alignment algorithm will
//' use optimizations that will cause it to exit early if the optimal alignment
//' has a distance greater than the distance threshold. This should not change
//' the correctness of distance calculations below the threshold, and results in
//' a large speedup. It is recommended to use `constrain=FALSE` only to verify
//' that the results do not change.
//' @param threads (`integer` count) number of parallel threads to use for
//' computation.
//' @param verbose (`integer` level) verbosity level
//'
//' @return (`data.frame`) a sparse distance matrix; columns are `seq1` and
//' `seq2` for the 0-based indices of two sequences; `score1` and `score2` are
//' the optimal alignment score for the two sequences in the "prealignment" (if
//' any) and "alignment" stages; `dist1` and `dist2` are the corresponding
//' sequence distances.
//'
//' @export
//' @rdname seq_distmx
// [[Rcpp::export]]
Rcpp::DataFrame seq_distmx_edlib(
  std::vector<std::string> seq,
  double dist_threshold,
  int details = 0,
  int span = 0,
  bool constrain = true,
  std::uint8_t threads = 1,
  int verbose = 0
) {

  std::vector<size_t> seq1, seq2;
  std::vector<int> score1, score2;
  std::vector<double> dist1, dist2;
  std::vector<int> aln_len, n_ins, n_del, max_ins, max_del;
  std::vector<std::string> cigar;

  AlignmentSpan span_enum;
  switch (span) {
    case 0: span_enum = AlignmentSpan::GLOBAL; break;
    case 1: span_enum = AlignmentSpan::EXTEND; break;
    default:
      Rcpp::stop("invalid value for 'span'");
  }

  std::unique_ptr<SparseDistanceMatrix> sdm;
  if (details == 0) {
    sdm = std::make_unique<SparseDistanceMatrix>(seq1, seq2, score1, score2, dist1, dist2);
  } else if (details == 1) {
    sdm = std::make_unique<SparseDistanceMatrixGapstats>(seq1, seq2, score1, score2, dist1, dist2, aln_len, n_ins, n_del, max_ins, max_del);
  } else if (details == 2) {
    sdm = std::make_unique<SparseDistanceMatrixCigar>(seq1, seq2, score1, score2, dist1, dist2, cigar);
  } else {
    Rcpp::stop("invalid value for 'details'");
  }
  std::unique_ptr<EdlibDistWorker> worker = create_edlib_dist_worker(seq,
                          dist_threshold, threads,
                          *sdm, verbose, span_enum, constrain);
  if (threads > 1) {
    RcppParallel::parallelFor(0, threads, *worker, 1, threads);
  } else {
    (*worker)(0, 1);
  }

  if (verbose >= 1) {
    Rcpp::Rcout << seq1.size() << " included / "
                << worker->aligned() << " aligned / "
                << worker->prealigned() << " prealigned"
                << std::endl;
  }

  Rcpp::DataFrame out = Rcpp::DataFrame::create(
    Rcpp::Named("seq1") = Rcpp::wrap(seq1),
    Rcpp::Named("seq2") = Rcpp::wrap(seq2),
    Rcpp::Named("score1") = Rcpp::wrap(score1),
    Rcpp::Named("score2") = Rcpp::wrap(score2),
    Rcpp::Named("dist1") = Rcpp::wrap(dist1),
    Rcpp::Named("dist2") = Rcpp::wrap(dist2)
  );
  return out;
}
