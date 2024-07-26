#include <cstdint>
#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>

extern "C" {
#include <SneakySnake.h>
}

#include "pairwise_alignment.h"
#include "pad_strings.h"
#include "SparseDistanceMatrix.h"

struct SneakySnakeAlignWorker : public RcppParallel::Worker {
  const std::vector<std::string> &seq;
  const std::shared_ptr<char[]> pseq;
  std::size_t seq_width = 0;
  const int match, mismatch, gap, extend, gap2, extend2;
  const double dist_threshold, sim_threshold, sim_threshold_plus_1;
  const bool is_constrained, is_score_constrained;
  const std::uint8_t threads;
  SparseDistanceMatrix &sdm;
  size_t &prealigned, &aligned;

  SneakySnakeAlignWorker(
    const std::vector<std::string> &seq,
    const int match,
    const int mismatch,
    const int gap,
    const int extend,
    const int gap2,
    const int extend2,
    const double dist_threshold,
    const bool constrain,
    const std::uint8_t threads,
    SparseDistanceMatrix &sdm,
    size_t &prealigned,
    size_t &aligned
  ) : seq(seq), pseq(pad_strings(seq, seq_width)),
  match(match), mismatch(mismatch), gap(gap), extend(extend),
  gap2(gap2), extend2(extend2), dist_threshold(dist_threshold),
  sim_threshold(1.0 - dist_threshold),
  sim_threshold_plus_1(1.0 + sim_threshold), is_constrained(constrain),
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

    wfa::WFAlignerChoose aligner{match, mismatch, gap, extend, gap2, extend2,
                                 wfa::WFAligner::Alignment};
    if (begin == 0) {
      begin_i = 1;
    } else {
      begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
    }
    size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
    sdm.mutex.lock();
    std::cout << "SneakySnakeAlignWorker thread " << begin << " entered; sequences [" <<
      begin_i << ", "<< end_i << ")" << std::endl;
    sdm.mutex.unlock();
    for (size_t i = begin_i; i < end_i; i++) {
      for (size_t j = 0; j < i; j++) {
        bool is_seqj_longer = seq[j].size() > seq[i].size();
        size_t s1 = is_seqj_longer ? i : j;
        size_t s2 = is_seqj_longer ? j : i;
        double l1 = seq[s1].size(), l2 = seq[s2].size();

        // sdm.mutex.lock();
        // Rcpp::Rcout << "#### seq " << i << " (l1=" << l1 << ") and "
        //             << j << " (l2=" << l2 <<")####" << std::endl;
        // sdm.mutex.unlock();

        if (l1/l2 < sim_threshold) continue;
        double maxd1 = dist_threshold * (l1 + l2) / sim_threshold_plus_1;

        int kmer_width = (int)ceil(2*(l2 - l1 * sim_threshold) / sim_threshold_plus_1);

        // sdm.mutex.lock();
        // Rcpp::Rcout << "maxd1=" << maxd1 << " kmer_width=" << kmer_width << std::endl;
        // sdm.mutex.unlock();

        int d1 = SneakySnake(
          l2,
          (char*)pseq.get() + s2*seq_width,
          (char*)pseq.get() + s1*seq_width,
          ceil(maxd1),
          kmer_width,
          0,
          l2
        );

        // sdm.mutex.lock();
        // Rcpp::Rcout << "SneakySnake: " << d1 << std::endl;
        // sdm.mutex.unlock();

        ++my_prealigned;
        if (d1 == 0) continue;
        // double d1 = 0;
        if (is_constrained) {
          int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
          int min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);
          aligner.setHeuristicBandedStatic(min_k, max_k);
          if (is_score_constrained) {
            aligner.setMaxAlignmentScore((int) maxd1 + 1);
          }
        }
        auto d2 = score_and_distance_wfa2(seq[s1], seq[s2], aligner);
        my_aligned++;
        if (d2.second <= dist_threshold) {
          my_seq1.push_back(j);
          my_seq2.push_back(i);
          my_score1.push_back(0);
          my_score2.push_back(d2.second);
          my_dist1.push_back(d1);
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

//' @param match (non-negative `integer`) alignment score for matching nucleotides
//' @param mismatch (non-negative `integer`) alignment penalty for mismatched
//' nucleotides.
//' @param gap_open (non-negative `integer`) alignment penalty for opening a new
//' gap (i.e., insertion or deletion).
//' @param gap_extend (non-negative `integer`) alignment penalty for each
//' position in a gap.
//' @param gap_open2 (non-negative `integer`) alternate alignment penalty for
//' opening a new gap (i.e., insertion or deletion).
//' @param gap_extend2 (non-negative `integer`) alternate alignment penalty for
//' each position in a gap.
//' @rdname seq_distmx
//' @export
// [[Rcpp::export]]
 Rcpp::DataFrame seq_distmx_snsn(
     const std::vector<std::string> &seq,
     const double dist_threshold,
     const int match = 1,
     const int mismatch = 2,
     const int gap_open = 10,
     const int gap_extend = 1,
     const int gap_open2 = 0,
     const int gap_extend2 = 0,
     const bool constrain = true,
     std::uint8_t threads = 1
 ) {
   size_t prealigned = 0, aligned = 0;

   std::vector<size_t> seq1, seq2;
   std::vector<int> score1, score2;
   std::vector<double> dist1, dist2;

   SparseDistanceMatrix sdm {seq1, seq2, score1, score2, dist1, dist2};
   SneakySnakeAlignWorker worker(seq,
                                 match, mismatch, gap_open, gap_extend, gap_open2,
                                 gap_extend2,
                                 dist_threshold, constrain, threads,
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
