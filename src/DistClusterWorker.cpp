#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>
#include <bindings/cpp/WFAligner.hpp>
#include <edlib.h>
#include "pairwise_alignment.h"
#include "ClusterAlgorithm.h"
#include "ClusterIndexedMatrix.h"

class HybridSplitClusterWorker : public RcppParallel::Worker {

  const std::vector<std::string> &seq;
  const double breakpoint;
  ClusterAlgorithm &clust_algo;
  const uint8_t threads;
  size_t &prealigned, &aligned;
  std::mutex mutex;
public :
  HybridSplitClusterWorker(
    const std::vector<std::string> &seq,
    const double breakpoint,
    ClusterAlgorithm &clust_algo,
    const uint8_t threads,
    size_t &prealigned,
    size_t &aligned
  ) : seq(seq), breakpoint(breakpoint), clust_algo(clust_algo),
  threads(threads), prealigned(prealigned), aligned(aligned) {};

  void operator()(std::size_t begin, std::size_t end) {
    double n = seq.size();
    double m = (n*n - 3.0*n + 2.0)/2.0;
    size_t my_prealigned = 0;
    size_t my_aligned = 0;
    size_t begin_i;

    EdlibAlignConfig ed_aligner = edlibNewAlignConfig(-1, EdlibAlignMode::EDLIB_MODE_NW, EdlibAlignTask::EDLIB_TASK_PATH, 0, 0);
    wfa::WFAlignerEdit wfa_aligner{wfa::WFAligner::Alignment};
    ClusterAlgorithm * my_algo = clust_algo.make_child();

    if (begin == 0) {
      begin_i = 1;
    } else {
      begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
    }
    size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
    mutex.lock();
    std::cout << "HybridSplit thread " << begin << " entered; sequences [" <<
      begin_i << ", "<< end_i << ")" << std::endl;
    mutex.unlock();
    for (size_t i = begin_i; i < end_i; i++) {
      for (size_t j = 0; j < i; j++) {
        double threshold = my_algo->max_relevant(i, j);
        // RcppThread::Rcout << "seqs " << j
        //                   << " and " << i
        //                   << " max relevant=" << threshold
        //                   << std::endl;

        bool is_seqj_longer = seq[j].size() > seq[i].size();
        size_t s1 = is_seqj_longer ? i : j;
        size_t s2 = is_seqj_longer ? j : i;
        double l1 = seq[s1].size(), l2 = seq[s2].size();
        // Rcpp::Rcout << "#### seq " << i << " (l1=" << l1 << ") and "
        //             << j << " (l2=" << l2 <<")####" << std::endl;

        double sim_threshold = 1.0 - threshold; // compiler can probably do this?
        if (l1/l2 < sim_threshold) continue;
        ++prealigned;
        double sim_threshold_plus_1 = 2.0 - threshold;
        double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
        bool is_close = breakpoint >= 1 ? maxd1 < breakpoint : threshold < breakpoint;
        double d;
        if (is_close) {
          int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
          int min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);
          wfa_aligner.setHeuristicBandedStatic(min_k, max_k);
          wfa_aligner.setMaxAlignmentScore((int)maxd1 + 1);
          d = distance_wfa2(seq[s1], seq[s2], wfa_aligner);
        } else {
          ed_aligner.k = (int)maxd1 + 1;
          d = distance_edlib(seq[s1], seq[s2], ed_aligner);
        }
        if (d < 1.0) ++my_aligned;

        // std::cout << "distance=" << d
        //                   << std::endl;
        if (d < threshold) (*my_algo)(j, i, d);
        RcppThread::checkUserInterrupt();
      }
    }
    mutex.lock();
    aligned += my_aligned;
    prealigned += my_prealigned;
    mutex.unlock();
    my_algo->merge_into_parent();
    std::cout << "thread " << begin << "done" << std::endl;
  }
};

Rcpp::IntegerMatrix single_linkage_hybrid_split(
    const std::vector<std::string> &seq,
    const DistanceConverter &dconv,
    const d_t m,
    const double breakpoint = 0.1,
    const uint8_t threads = 1
) {
  size_t prealigned = 0, aligned = 0;

  Rcpp::IntegerMatrix im(m, seq.size());
  ClusterIndexedMatrix<RcppParallel::RMatrix<int>> cm(dconv, im);
  HybridSplitClusterWorker worker(seq, breakpoint, cm, threads, prealigned, aligned);
  if (threads > 1) {
    RcppParallel::parallelFor(0, threads, worker, 1, threads);
  } else {
    worker(0, 1);
  }
  std::cout << aligned << " aligned / "
            << prealigned << " prealigned"
            << std::endl;
  return im;
}

//' @export
 // [[Rcpp::export]]
 Rcpp::IntegerMatrix single_linkage_hybrid_split_uniform(
     const std::vector<std::string> &seq,
     const float dmin,
     const float dmax,
     const float dstep,
     const double breakpoint = 0.1,
     const int threads = 1
 ) {
   const UniformDistanceConverter dconv(dmin, dmax, dstep);
   const int m = (int) ceilf((dmax - dmin)/dstep) + 1;
   return single_linkage_hybrid_split(seq, dconv, m, breakpoint, threads);
 }

//' @export
 // [[Rcpp::export]]
 Rcpp::IntegerMatrix single_linkage_hybrid_split_array(
     const std::vector<std::string> &seq,
     const std::vector<double> &thresholds,
     const double breakpoint = 0.1,
     const int threads = 1
 ) {
   const ArrayDistanceConverter dconv(thresholds);
   const int m = thresholds.size();
   return single_linkage_hybrid_split(seq, dconv, m, breakpoint, threads);
 }

//' @export
 // [[Rcpp::export]]
 Rcpp::IntegerMatrix single_linkage_hybrid_split_cached(
     const std::vector<std::string> &seq,
     const std::vector<double> &thresholds,
     const double precision,
     const double breakpoint = 0.1,
     const int threads = 1
 ) {
   const CachedDistanceConverter dconv(thresholds, precision);
   const int m = thresholds.size();
   return single_linkage_hybrid_split(seq, dconv, m, breakpoint, threads);
 }

class HybridConcurrentClusterWorker : public RcppParallel::Worker {

  const std::vector<std::string> &seq;
  const double breakpoint;
  ClusterAlgorithm &clust_algo;
  const uint8_t threads;
  size_t &prealigned, &aligned;
  std::mutex mutex;
public :
  HybridConcurrentClusterWorker(
    const std::vector<std::string> &seq,
    const double breakpoint,
    ClusterAlgorithm &clust_algo,
    const uint8_t threads,
    size_t &prealigned,
    size_t &aligned
  ) : seq(seq), breakpoint(breakpoint), clust_algo(clust_algo),
  threads(threads), prealigned(prealigned), aligned(aligned) {};

  void operator()(std::size_t begin, std::size_t end) {
    double n = seq.size();
    double m = (n*n - 3.0*n + 2.0)/2.0;
    size_t my_prealigned = 0;
    size_t my_aligned = 0;
    size_t begin_i;

    EdlibAlignConfig ed_aligner = edlibNewAlignConfig(
      -1,
      EdlibAlignMode::EDLIB_MODE_NW,
      EdlibAlignTask::EDLIB_TASK_PATH,
      0,
      0
    );
    wfa::WFAlignerEdit wfa_aligner{wfa::WFAligner::Alignment};

    if (begin == 0) {
      begin_i = 1;
    } else {
      begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
    }
    size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
    mutex.lock();
    std::cout << "HybridConcurrent thread " << begin << " entered; sequences [" <<
      begin_i << ", "<< end_i << ")" << std::endl;
    mutex.unlock();
    for (size_t i = begin_i; i < end_i; i++) {
      for (size_t j = 0; j < i; j++) {
        // mutex.lock();
        // std::cout << "Thread " << begin
        //           << ": seqs " << j
        //           << " and " << i
        //           << std::endl;
        // mutex.unlock();
        double threshold = clust_algo.max_relevant(i, j);
        // mutex.lock();
        // std::cout << "Thread " << begin
        //           << ": max relevant=" << threshold
        //                   << std::endl;
        // mutex.unlock();
        bool is_seqj_longer = seq[j].size() > seq[i].size();
        size_t s1 = is_seqj_longer ? i : j;
        size_t s2 = is_seqj_longer ? j : i;
        double l1 = seq[s1].size(), l2 = seq[s2].size();
        // Rcpp::Rcout << "#### seq " << i << " (l1=" << l1 << ") and "
        //             << j << " (l2=" << l2 <<")####" << std::endl;

        double sim_threshold = 1.0 - threshold; // compiler can probably do this?
        if (l1/l2 < sim_threshold) continue;
        ++my_prealigned;
        double sim_threshold_plus_1 = 2.0 - threshold;
        double maxd1 = threshold * (l1 + l2) / sim_threshold_plus_1;
        bool is_close = breakpoint >= 1 ? maxd1 < breakpoint : threshold < breakpoint;
        double d;
        if (is_close) {
          int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
          int min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);
          wfa_aligner.setHeuristicBandedStatic(min_k, max_k);
          wfa_aligner.setMaxAlignmentScore((int)maxd1 + 1);
          // std::cout << "wfa_aligner min_k=" << min_k
          //           << " max_k=" << max_k
          //           << " max score=" << (int)maxd1 + 1
          //           << std::endl;
          d = distance_wfa2(seq[s1], seq[s2], wfa_aligner);
        } else {
          ed_aligner.k = (int)maxd1 + 1;
          d = distance_edlib(seq[s1], seq[s2], ed_aligner);
        }
        if (d < 1.0) ++my_aligned;
        // mutex.lock();
        // std::cout << "Thread " << begin
        //           << ": distance=" << d
        //           << std::endl;
        // mutex.unlock();
        if (d < threshold) clust_algo(j, i, d);
        mutex.lock();
        // std::cout << "Thread " << begin
        //           << ": finished " << j
        //           << " and " << i
        //           << std::endl;
        mutex.unlock();
        RcppThread::checkUserInterrupt();
      }
    }
    mutex.lock();
    aligned += my_aligned;
    prealigned += my_prealigned;
    std::cout << "Exiting thread " << begin << std::endl;
    mutex.unlock();
  }
};

Rcpp::IntegerMatrix single_linkage_hybrid_concurrent(
    const std::vector<std::string> &seq,
    const DistanceConverter &dconv,
    const d_t m,
    const double breakpoint = 0.1,
    const uint8_t threads = 1
) {
  size_t prealigned = 0, aligned = 0;

  Rcpp::IntegerMatrix im(m, seq.size());
  ClusterIndexedMatrix<RcppParallel::RMatrix<int>> cm(dconv, im);
  HybridConcurrentClusterWorker worker(seq, breakpoint, cm, threads, prealigned, aligned);
  if (threads > 1) {
    RcppParallel::parallelFor(0, threads, worker, 1, threads);
  } else {
    worker(0, 1);
  }
  std::cout << "back in main function" << std::endl;
  std::cout << aligned << " aligned / "
            << prealigned << " prealigned"
            << std::endl;
  return im;
}

//' @export
 // [[Rcpp::export]]
 Rcpp::IntegerMatrix single_linkage_hybrid_concurrent_uniform(
     const std::vector<std::string> &seq,
     const float dmin,
     const float dmax,
     const float dstep,
     const double breakpoint = 0.1,
     const int threads = 1
 ) {
   const UniformDistanceConverter dconv(dmin, dmax, dstep);
   const int m = (int) ceilf((dmax - dmin)/dstep) + 1;
   return single_linkage_hybrid_concurrent(seq, dconv, m, breakpoint, threads);
 }

//' @export
 // [[Rcpp::export]]
 Rcpp::IntegerMatrix single_linkage_hybrid_concurrent_array(
     const std::vector<std::string> &seq,
     const std::vector<double> &thresholds,
     const double breakpoint = 0.1,
     const int threads = 1
 ) {
   const ArrayDistanceConverter dconv(thresholds);
   const int m = thresholds.size();
   return single_linkage_hybrid_concurrent(seq, dconv, m, breakpoint, threads);
 }

//' @export
 // [[Rcpp::export]]
 Rcpp::IntegerMatrix single_linkage_hybrid_concurrent_cached(
     const std::vector<std::string> &seq,
     const std::vector<double> &thresholds,
     const double precision,
     const double breakpoint = 0.1,
     const int threads = 1
 ) {
   const CachedDistanceConverter dconv(thresholds, precision);
   const int m = thresholds.size();
   return single_linkage_hybrid_concurrent(seq, dconv, m, breakpoint, threads);
 }
