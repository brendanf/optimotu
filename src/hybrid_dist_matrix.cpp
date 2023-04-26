#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>

extern "C" {
#include <SneakySnake.h>
}

#include "pairwise_alignment.h"
#include "pad_strings.h"
#include "SparseDistanceMatrix.h"

struct HybridAlignWorker : public RcppParallel::Worker {
  const std::vector<std::string> &seq;
  const std::shared_ptr<char[]> pseq;
  const double dist_threshold, sim_threshold, sim_threshold_plus_1;
  const double breakpoint;
  const uint8_t threads;
  SparseDistanceMatrix &sdm;
  size_t &prealigned, &aligned;
  size_t seq_width = 0;

  HybridAlignWorker(
    const std::vector<std::string> &seq,
    const double dist_threshold,
    const double breakpoint,
    const uint8_t threads,
    SparseDistanceMatrix &sdm,
    size_t &prealigned,
    size_t &aligned
  ) : seq(seq), pseq(pad_strings(seq, seq_width)),
  dist_threshold(dist_threshold), sim_threshold(1.0 - dist_threshold),
  sim_threshold_plus_1(1.0 + sim_threshold),
  breakpoint(breakpoint),
  threads(threads),
  sdm(sdm), prealigned(prealigned), aligned(aligned) {};

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

    EdlibAlignConfig ed_aligner = edlibNewAlignConfig(-1, EdlibAlignMode::EDLIB_MODE_NW, EdlibAlignTask::EDLIB_TASK_PATH, 0, 0);
    wfa::WFAlignerEdit wfa_aligner{wfa::WFAligner::Alignment};

    if (begin == 0) {
      begin_i = 1;
    } else {
      begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
    }
    size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
    sdm.mutex.lock();
    std::cout << "HybridAlignWorker thread " << begin << " entered; sequences [" <<
      begin_i << ", "<< end_i << ")" << std::endl;
    sdm.mutex.unlock();
    for (size_t i = begin_i; i < end_i; i++) {
      for (size_t j = 0; j < i; j++) {

        bool is_seqj_longer = seq[j].size() > seq[i].size();
        size_t s1 = is_seqj_longer ? i : j;
        size_t s2 = is_seqj_longer ? j : i;
        double l1 = seq[s1].size(), l2 = seq[s2].size();
        // Rcpp::Rcout << "#### seq " << i << " (l1=" << l1 << ") and "
        //             << j << " (l2=" << l2 <<")####" << std::endl;

        if (l1/l2 < sim_threshold) continue;
        double maxd1 = dist_threshold * (l1 + l2) / sim_threshold_plus_1;
        bool is_close = false;
        std::pair<int, double> d1 = {0, 0};
        if (breakpoint >= 1) {
          // edit distance, not pairwise dissimilarity
          is_close = SneakySnake(
            l2,
            pseq.get() + s2*seq_width,
            pseq.get() + s1*seq_width,
            (int) breakpoint,
            (int) breakpoint * 2 + 1,
            0,
            l2
          );
          ++my_prealigned;
          if (!is_close && maxd1 <= breakpoint) continue;
        } else {
          int kmer_width = (int)round(2*(l2 - l1 * (1-breakpoint)) /
                            (2-breakpoint));
          is_close = SneakySnake(
            l2,
            pseq.get() + s2*seq_width,
            pseq.get() + s1*seq_width,
            (int) (breakpoint * (l1 + l2) / (2 - breakpoint)),
            kmer_width,
            0,
            l2
          );
          ++my_prealigned;
          if (!is_close && dist_threshold <= breakpoint) continue;
        }
        int sc2;
        double d2;
        if (is_close) {
          int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
          int min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);
          wfa_aligner.setHeuristicBandedStatic(min_k, max_k);
          wfa_aligner.setMaxAlignmentScore((int)maxd1 + 1);
          auto sd2 = score_and_distance_wfa2(seq[s1], seq[s2], wfa_aligner);
          sc2 = sd2.first;
          d2 = sd2.second;
          my_aligned++;
        } else {
          ed_aligner.k = (int)maxd1 + 1;

          auto alignResult = edlibAlign(
            seq[s1].c_str(), seq[s1].size(),
            seq[s2].c_str(), seq[s2].size(),
            ed_aligner
          );
          my_aligned++;
          if (alignResult.status == EDLIB_STATUS_ERROR) continue;
          if (alignResult.editDistance == -1) continue;
          sc2 = alignResult.editDistance;
          d2 = (double)alignResult.editDistance / (double) alignResult.alignmentLength;
          edlibFreeAlignResult(alignResult);
        }
        if (d2 <= dist_threshold) {
          my_seq1.push_back(j);
          my_seq2.push_back(i);
          my_score1.push_back(d1.first);
          my_score2.push_back(sc2);
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

//' @param breakpoint (`numeric` scalar) threshold for deciding whether to use
//' WFA2 or edlib for edit-distance alignment.  This parameter is interpreted as
//' an edit distance if greater than or equal to `1`, or as a pairwise
//' dissimilarity if less than 1. In either case, WFA2 is used below the
//' breakpoint, and edlib is used above it.
//' @export
//' @rdname seq_distmx
// [[Rcpp::export]]
Rcpp::DataFrame seq_distmx_hybrid(std::vector<std::string> seq, double dist_threshold,
                                  double breakpoint = 0.1, uint8_t threads = 1) {
  size_t prealigned = 0, aligned = 0;

  std::vector<size_t> seq1, seq2;
  std::vector<int> score1, score2;
  std::vector<double> dist1, dist2;

  SparseDistanceMatrix sdm {seq1, seq2, score1, score2, dist1, dist2};
  HybridAlignWorker worker(seq,
                           dist_threshold, breakpoint, threads,
                           sdm, prealigned, aligned);
  if (threads > 1) {
    RcppParallel::parallelFor(0, threads, worker, 1, threads);
  } else {
    worker(0, 1);
  }

  Rcpp::Rcout << seq1.size() << " included / "
              << aligned << " aligned / "
              << prealigned << " prealigned"
              << std::endl;

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
