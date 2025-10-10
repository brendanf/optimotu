#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>

#include "SparseDistanceMatrix.h"
#include "PackedSequenceSet.h"

struct HammingAlignWorker : public RcppParallel::Worker {
  const PackedSequenceSet seq;
  const double dist_threshold;
  const int min_overlap;
  const bool ignore_gap;
  const std::uint8_t threads;
  SparseDistanceMatrix &sdm;
  size_t &aligned;
  size_t &considered;
  const int verbose;

  HammingAlignWorker(
    const std::vector<std::string> &seq,
    const double dist_threshold,
    const int min_overlap,
    const bool ignore_gap,
    const std::uint8_t threads,
    SparseDistanceMatrix &sdm,
    size_t &aligned,
    size_t &considered,
    const int verbose
  ) :
    seq(seq),
    dist_threshold(dist_threshold),
    min_overlap(min_overlap),
    ignore_gap(ignore_gap),
    threads(threads),
    sdm(sdm),
    aligned(aligned),
    considered(considered),
    verbose(verbose) {};

  void operator()(std::size_t begin, std::size_t end) {
    double n = seq.num_seqs;
    double m = (n*n - 3.0*n + 2.0)/2.0;
    size_t my_aligned = 0;
    size_t my_considered = 0;
    size_t begin_i;
    if (begin == 0) {
      begin_i = 1;
    } else {
      begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
    }
    size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
    if (verbose > 1) {
      sdm.mutex.lock();
      RcppThread::Rcout << "HammingAlignWorker thread " << begin
                        << " entered; sequences [" << begin_i
                        << ", "<< end_i << ")" << std::endl;
      sdm.mutex.unlock();
    }
    std::vector<size_t> my_seq1;
    my_seq1.reserve(1000);
    std::vector<size_t> my_seq2;
    my_seq2.reserve(1000);
    std::vector<int> my_score1;
    my_score1.reserve(1000);
    std::vector<int> my_score2;
    std::vector<double> my_dist1;
    my_dist1.reserve(1000);
    std::vector<double> my_dist2;
    my_dist2.reserve(1000);

    for (size_t i = begin_i; i < end_i; i++) {
      for (size_t j = 0; j < i; j++) {
        if (verbose > 2) {
          Rcpp::Rcerr << "considering " << i
                      << " (seq: " << seq.packed_seq[i].size()
                      << ", mask: " << seq.mask[i].size()
                      << ") and " << j
                      << " (seq: " << seq.packed_seq[j].size()
                      << ", mask: " << seq.mask[j].size()
                      << ")"<< std::endl;
        }
        auto[success, num_matches, d] =
          seq.success_score_and_dist(i, j, min_overlap, ignore_gap);
        my_considered++;
        if (success) {
          my_aligned++;
          if (d <= dist_threshold) {
            my_seq1.push_back(j);
            my_seq2.push_back(i);
            my_score1.push_back(0);
            my_score2.push_back(num_matches);
            my_dist1.push_back(0);
            my_dist2.push_back(d);
            if (my_seq1.size() == 1000) {
              sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2);
            }
          }
        }
      }
    }
    if (my_seq1.size() > 0) {
      sdm.append(my_seq1, my_seq2, my_score1, my_score2, my_dist1, my_dist2);
    }
    sdm.mutex.lock();
    aligned += my_aligned;
    considered += my_considered;
    sdm.mutex.unlock();
  }
};

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
  size_t aligned = 0;
  size_t considered = 0;
  HammingAlignWorker worker(seq, dist_threshold, min_overlap, ignore_gap,
    threads, sdm, aligned, considered, verbose);
  worker.seq.verify(seq);
  if (verbose >= 2) {
    Rcpp::Rcerr << "finished initializing worker" << std::endl;
  }
  if (threads == 1) {
    worker(0, 1);
  } else {
    RcppParallel::parallelFor(0, threads, worker);
  }
  if (verbose >= 1) {
    Rcpp::Rcerr << seq1.size() << " included / "
                << aligned << " aligned / "
                << considered << " considered" << std::endl;
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
