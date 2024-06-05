#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>

#include "pairwise_alignment.h"
#include "SparseDistanceMatrix.h"
#include "ClusterIndexedMatrix.h"
#include <cstdint>

struct PrealignAlignWorker : public RcppParallel::Worker {
  const std::vector<std::string> &seq;
  const int match, mismatch, gap, extend, gap2, extend2;
  const double dist_threshold, sim_threshold, sim_threshold_plus_1;
  const bool do_prealign, is_constrained, is_score_constrained;
  const std::uint8_t threads;
  SparseDistanceMatrix &sdm;
  size_t &prealigned, &aligned;

  PrealignAlignWorker(
    const std::vector<std::string> &seq,
    const int match,
    const int mismatch,
    const int gap,
    const int extend,
    const int gap2,
    const int extend2,
    const double dist_threshold,
    const bool do_prealign,
    const bool constrain,
    const std::uint8_t threads,
    SparseDistanceMatrix &sdm,
    size_t &prealigned,
    size_t &aligned
  ) : seq(seq),
  match(match), mismatch(mismatch), gap(gap), extend(extend), gap2(gap2), extend2(extend2),
  dist_threshold(dist_threshold), sim_threshold(1.0 - dist_threshold),
  sim_threshold_plus_1(1.0 + sim_threshold), do_prealign(do_prealign),
  is_constrained(constrain),
  is_score_constrained(
    constrain &&
      match == 0 &&
      mismatch == 1 &&
      gap == 0 &&
      extend == 1 &&
      gap2 == 0 &&
      (extend2 == 0 || extend2 == 1)
  ),
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

    wfa::WFAlignerEdit prealigner{wfa::WFAligner::Score};

    wfa::WFAlignerChoose aligner{match, mismatch, gap, extend, gap2, extend2,
                                 wfa::WFAligner::Alignment};
    if (begin == 0) {
      begin_i = 1;
    } else {
      begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
    }
    size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
    sdm.mutex.lock();
    std::cout << "PrealignAlignWorker thread " << begin << " entered; sequences [" <<
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
        int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
        int min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);

        std::pair<int, double> d1 = {0, 0};
        if (do_prealign) {
          // std::cout << "Setting max score to " << min_k << ", " << max_k << std::endl;
          prealigner.setMaxAlignmentScore((int) maxd1 + 1);
          // std::cout << "Setting band heuristics to " << min_k << ", " << max_k << std::endl;
          prealigner.setHeuristicBandedStatic(min_k, max_k);
          // std::cout << "Prealigning..." << std::endl;
          auto status = prealigner.alignEnd2End(seq[s1], seq[s2]);
          // std::cout << "Prealignment finished." << std::endl;
          ++my_prealigned;
          if (status != wfa::WFAligner::AlignmentStatus::StatusSuccessful) continue;
          // std::cout << "Prealignment successful." << std::endl;
          int sc1 = prealigner.getAlignmentScore();
          if (sc1 > maxd1) continue;
          d1 = {sc1, (double)sc1*sim_threshold_plus_1 / (l1 + l2)};
        }
        if (is_constrained) {
          // std::cout << "Setting band heuristics to " << min_k << ", " << max_k << std::endl;
          aligner.setHeuristicBandedStatic(min_k, max_k);
          if (is_score_constrained) {
            // std::cout << "Setting max score to " << (int)maxd1 + 1 << std::endl;
            aligner.setMaxAlignmentScore((int)maxd1 + 1);
          }
        }
        // std::cout << "Aligning..." << std::endl;
        auto d2 = score_and_distance_wfa2(seq[s1], seq[s2], aligner);
        // std::cout << "Alignment successful." << std::endl;

        my_aligned++;
        if (d2.second <= dist_threshold) {
          my_seq1.push_back(j);
          my_seq2.push_back(i);
          my_score1.push_back(d1.first);
          my_score2.push_back(d2.first);
          my_dist1.push_back(d1.second);
          my_dist2.push_back(d2.second);

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

//' @param prealign (`logical` flag) if `TRUE`, do a prealignment using
//' edit-distance as alignment score (`match = 0, mismatch = 1, gap_open = 0,
//' gap_extend = 1, gap_open2 = 0, gap_extend2 = 1`) to test feasibility before
//' aligning with alternate alignment scores. Note that, since pairwise distance
//' is defined using edit distance, any other set of scores will always result
//' in a pairwise distance which is equal to or greater than an alignment based
//' on the edit distance score.
//' @export
//' @rdname seq_distmx
// [[Rcpp::export]]
Rcpp::DataFrame seq_distmx_wfa2(
    std::vector<std::string> seq,
    double dist_threshold,
    int match = 1, int mismatch = 2,
    int gap_open = 10, int gap_extend = 1,
    int gap_open2 = 0, int gap_extend2 = 0,
    bool prealign = true,
    bool constrain = true,
    std::uint8_t threads = 1) {
  size_t prealigned = 0, aligned = 0;

  std::vector<size_t> seq1, seq2;
  std::vector<int> score1, score2;
  std::vector<double> dist1, dist2;

  SparseDistanceMatrix sdm {seq1, seq2, score1, score2, dist1, dist2};
  PrealignAlignWorker worker(seq,
                             match, mismatch, gap_open, gap_extend, gap_open2,
                             gap_extend2,
                             dist_threshold, prealign, constrain, threads,
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
