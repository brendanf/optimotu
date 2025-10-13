#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>

#include "SparseDistanceMatrix.h"
#include "PackedSequenceSet.h"

#include "HammingDistWorker.h"

//' @export
//' @rdname seq_distmx
// [[Rcpp::export]]
Rcpp::DataFrame seq_distmx_hamming(
  const std::vector<std::string> seq,
  const double dist_threshold,
  const int min_overlap = 0,
  const bool ignore_gap = false,
  const std::uint8_t threads = 1,
  const int verbose = 0
) {
  if (verbose >= 2) {
    Rcpp::Rcerr << "seq_distmx_hamming: " << seq.size() << " sequences" << std::endl;
  }
  std::vector<size_t> seq1, seq2;
  std::vector<int> score1, score2;
  std::vector<double> dist1, dist2;
  SparseDistanceMatrix sdm {seq1, seq2, score1, score2, dist1, dist2};
  std::unique_ptr<HammingDistWorker> worker = create_hamming_dist_worker(seq, dist_threshold,
    threads, sdm, min_overlap, ignore_gap, verbose);
  if (verbose >= 2) {
    Rcpp::Rcerr << "finished initializing worker" << std::endl;
  }
  if (threads == 1) {
    (*worker)(0, 1);
  } else {
    RcppParallel::parallelFor(0, threads, *worker);
  }
  if (verbose >= 1) {
    Rcpp::Rcerr << seq1.size() << " included / "
                << worker->aligned() << " aligned / "
                << worker->prealigned() << " considered" << std::endl;
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("seq1") = seq1,
    Rcpp::Named("seq2") = seq2,
    Rcpp::Named("score1") = score1,
    Rcpp::Named("score2") = score2,
    Rcpp::Named("dist1") = dist1,
    Rcpp::Named("dist2") = dist2
  );
}
