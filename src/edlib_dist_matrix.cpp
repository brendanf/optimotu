#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>

#include "pairwise_alignment.h"
#include "SparseDistanceMatrix.h"

#include <cstdint>

struct EdlibAlignWorker : public RcppParallel::Worker {
  const std::vector<std::string> &seq;
  const double dist_threshold, sim_threshold, sim_threshold_plus_1;
  const bool is_constrained;
  const std::uint8_t threads;
  SparseDistanceMatrix &sdm;
  size_t &prealigned, &aligned;
  const int verbose = 0;

  EdlibAlignWorker(
    const std::vector<std::string> &seq,
    const double dist_threshold,
    const bool constrain,
    const std::uint8_t threads,
    SparseDistanceMatrix &sdm,
    size_t &prealigned,
    size_t &aligned,
    const int verbose = 0
  ) : seq(seq),
  dist_threshold(dist_threshold), sim_threshold(1.0 - dist_threshold),
  sim_threshold_plus_1(1.0 + sim_threshold),
  is_constrained(constrain),
  threads(threads),
  sdm(sdm), prealigned(prealigned), aligned(aligned), verbose(verbose) {};

  void operator()(std::size_t begin, std::size_t end) {
    double n = seq.size();
    double m = (n*n - 3.0*n + 2.0)/2.0;
    size_t my_prealigned = 0;
    size_t my_aligned = 0;
    size_t begin_i;
    std::vector<size_t> my_seq1;
    my_seq1.reserve(1000);
    std::vector<size_t> my_seq2;
    my_seq2.reserve(1000);
    std::vector<int> my_score1;
    my_score1.reserve(1000);
    std::vector<int> my_score2;
    my_score2.reserve(1000);
    std::vector<double> my_dist1;
    my_dist1.reserve(1000);
    std::vector<double> my_dist2;
    my_dist2.reserve(1000);

    EdlibAlignConfig aligner = edlibNewAlignConfig(-1, EdlibAlignMode::EDLIB_MODE_NW, EdlibAlignTask::EDLIB_TASK_PATH, 0, 0);

    if (begin == 0) {
      begin_i = 1;
    } else {
      begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
    }
    size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));

    if (verbose >= 2) {
      RcppThread::Rcerr << "EdlibAlignWorker thread " << begin << " entered; sequences [" <<
        begin_i << ", "<< end_i << ")" << std::endl;
    }
    for (size_t i = begin_i; i < end_i; i++) {
      for (size_t j = 0; j < i; j++) {

        bool is_seqj_longer = seq[j].size() > seq[i].size();
        size_t s1 = is_seqj_longer ? i : j;
        size_t s2 = is_seqj_longer ? j : i;
        double l1 = seq[s1].size(), l2 = seq[s2].size();
        if (verbose >= 3) {
          RcppThread::Rcerr << "#### seq " << i << " (l1=" << l1 << ") and "
                            << j << " (l2=" << l2 <<")####" << std::endl;
        }

        if (l1/l2 < sim_threshold) continue;

        std::pair<int, double> d1 = {0, 0};
        if (is_constrained) {
          double maxd1 = dist_threshold * (l1 + l2) / sim_threshold_plus_1;
          aligner.k = (int)maxd1 + 1;
        }
        auto alignResult = edlibAlign(
          seq[s1].c_str(), seq[s1].size(),
          seq[s2].c_str(), seq[s2].size(),
          aligner
        );
        my_prealigned++;
        if (alignResult.status == EDLIB_STATUS_ERROR) continue;
        if (alignResult.editDistance == -1) continue;
        my_aligned++;
        double d2 = (double)alignResult.editDistance / (double) alignResult.alignmentLength;
        edlibFreeAlignResult(alignResult);
        if (d2 <= dist_threshold) {
          my_seq1.push_back(j);
          my_seq2.push_back(i);
          my_score1.push_back(d1.first);
          my_score2.push_back(alignResult.editDistance);
          my_dist1.push_back(d1.second);
          my_dist2.push_back(d2);

          if (my_seq1.size() == 1000) {
            sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2);
          }
        }
        RcppThread::checkUserInterrupt();
      }
    }
    if (my_seq1.size() > 0) {
      sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2);
    }
    sdm.mutex.lock();
    aligned += my_aligned;
    prealigned += my_prealigned;
    sdm.mutex.unlock();
  }
};

//' Sparse distance matrix between DNA sequences
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
Rcpp::DataFrame seq_distmx_edlib(std::vector<std::string> seq, double dist_threshold,
                                 bool constrain = true, std::uint8_t threads = 1,
                                 int verbose = 0) {
  size_t prealigned = 0, aligned = 0;

  std::vector<size_t> seq1, seq2;
  std::vector<int> score1, score2;
  std::vector<double> dist1, dist2;

  SparseDistanceMatrix sdm {seq1, seq2, score1, score2, dist1, dist2};
  EdlibAlignWorker worker(seq,
                          dist_threshold, constrain, threads,
                          sdm, prealigned, aligned, verbose);
  if (threads > 1) {
    RcppParallel::parallelFor(0, threads, worker, 1, threads);
  } else {
    worker(0, 1);
  }

  if (verbose >= 1) {
    Rcpp::Rcout << seq1.size() << " included / "
                << aligned << " aligned / "
                << prealigned << " prealigned"
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
