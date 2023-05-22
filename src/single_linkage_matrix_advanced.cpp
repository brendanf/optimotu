#ifdef OPTIMOTU_R
#include <fstream>
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "ClusterMatrix.h"
#include "ClusterIndexedMatrix.h"

std::unique_ptr<ClusterAlgorithm> create_clustermatrix(
    const DistanceConverter &dconv,
    Rcpp::IntegerMatrix im,
    const bool do_binary_search = true,
    const bool do_binary_fill = true,
    const bool do_topdown_fill = false
) {
  if (do_binary_search) {
    if (do_binary_fill) {
      return std::make_unique<ClusterMatrix<RcppParallel::RMatrix<int>, true, true, false>>(dconv, im);
    } else if (do_topdown_fill) {
      return std::make_unique<ClusterMatrix<RcppParallel::RMatrix<int>, true, false, true>>(dconv, im);
    } else {
      return std::make_unique<ClusterMatrix<RcppParallel::RMatrix<int>, true, false, false>>(dconv, im);
    }
  } else {
    if (do_binary_fill) {
      return std::make_unique<ClusterMatrix<RcppParallel::RMatrix<int>, false, true, false>>(dconv, im);
    } else if (do_topdown_fill) {
      return std::make_unique<ClusterMatrix<RcppParallel::RMatrix<int>, false, false, true>>(dconv, im);
    } else {
      return std::make_unique<ClusterMatrix<RcppParallel::RMatrix<int>, false, false, false>>(dconv, im);
    }
  }
}

class ClusterWorker : public RcppParallel::Worker {
protected:
  std::vector<ClusterAlgorithm*> algo_list;
  std::istream &file;
  std::mutex mutex;
public:
  ClusterWorker(ClusterAlgorithm *algo, std::istream &file, const int threads) :
    algo_list(threads),
    file(file) {
    std::cout << "ClusterWorker constructor start...";
    for (int i = 0; i < threads; ++i) {
      algo_list[i] = algo->make_child();
    }
    std::cout << "done" << std::endl;
  };

  void operator()(size_t begin, size_t end) {
    DistanceElement d;
    std::vector<DistanceElement> buffer;
    buffer.reserve(100);
    std::cout << "Starting ClusterWorker thread " << begin << std::endl;
    while (true) {
      mutex.lock();
      for (int i = 0; i < 100 && file; ++i) {
        file >> d;
        buffer.push_back(d);
      }
      mutex.unlock();
      if (buffer.size() == 0) break;
      for (auto d : buffer) {
        algo_list[begin]->operator()(d);
      }
      buffer.clear();
    }
    mutex.lock();
    std::cout << "ClusterWorker thread " << begin << " merging..." << std::endl;
    mutex.unlock();
    algo_list[begin]->merge_into_parent();

    mutex.lock();
    std::cout << "ClusterWorker thread " << begin << " exiting" << std::endl;
    mutex.unlock();
  }
};

Rcpp::IntegerMatrix single_linkage_matrix2(
    const std::string file,
    const Rcpp::CharacterVector &seqnames,
    const DistanceConverter &dconv,
    const int m,
    const int threads = 1,
    const bool do_binary_search = true,
    const bool do_binary_fill = true,
    const bool do_topdown_fill = false
) {
  Rcpp::IntegerMatrix im(m, seqnames.size());
  auto cm = create_clustermatrix(dconv, im, do_binary_search, do_binary_fill, do_topdown_fill);
  std::ifstream infile(file);
  if (!infile.is_open()) {
    Rcpp::stop("failed to open input file");
  }
  if (threads == 1) {
    DistanceElement d;
    while(infile >> d) {
      cm->operator()(d);
    }
  } else {
    ClusterWorker cw(cm.get(), infile, threads);
    RcppParallel::parallelFor(0, threads, cw, 1, threads);
  }
  infile.close();
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
      const int threads = 1,
      const bool do_binary_search = true,
      const bool do_binary_fill = true,
      const bool do_topdown_fill = false
) {
   const UniformDistanceConverter dconv(dmin, dstep);
   const int m = (int) ceilf((dmax - dmin)/dstep) + 1;
   return single_linkage_matrix2(file, seqnames, dconv, m, threads,
                                 do_binary_search,
                                 do_binary_fill, do_topdown_fill);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix single_linkage_matrix2_array(
      const std::string file,
      const Rcpp::CharacterVector &seqnames,
      const std::vector<double> &thresholds,
      const int threads = 1,
      const bool do_binary_search = true,
      const bool do_binary_fill = true,
      const bool do_topdown_fill = false
) {
   const ArrayDistanceConverter dconv(thresholds);
   const int m = thresholds.size();
   return single_linkage_matrix2(file, seqnames, dconv, m, threads,
                                 do_binary_search, do_binary_fill, do_topdown_fill);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix single_linkage_matrix2_cached(
      const std::string file,
      const Rcpp::CharacterVector &seqnames,
      const std::vector<double> &thresholds,
      const double precision,
      const int threads = 1,
      const bool do_binary_search = true,
      const bool do_binary_fill = true,
      const bool do_topdown_fill = false
) {
   const CachedDistanceConverter dconv(thresholds, precision);
   const int m = thresholds.size();
   return single_linkage_matrix2(file, seqnames, dconv, m, threads,
                                 do_binary_search,
                                 do_binary_fill, do_topdown_fill);
}

Rcpp::IntegerMatrix single_linkage_imatrix(
    const std::string file,
    const Rcpp::CharacterVector &seqnames,
    const DistanceConverter &dconv,
    const int m
) {
  Rcpp::IntegerMatrix im(m, seqnames.size());
  ClusterIndexedMatrix<RcppParallel::RMatrix<int>> cm(dconv, im);
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
 Rcpp::IntegerMatrix single_linkage_imatrix_uniform(
     const std::string file,
     const Rcpp::CharacterVector &seqnames,
     const float dmin,
     const float dmax,
     const float dstep
 ) {
   const UniformDistanceConverter dconv(dmin, dstep);
   const int m = (int) ceilf((dmax - dmin)/dstep) + 1;
   return single_linkage_imatrix(file, seqnames, dconv, m);
 }

//' @export
 // [[Rcpp::export]]
 Rcpp::IntegerMatrix single_linkage_imatrix_array(
     const std::string file,
     const Rcpp::CharacterVector &seqnames,
     const std::vector<double> &thresholds
 ) {
   const ArrayDistanceConverter dconv(thresholds);
   const int m = thresholds.size();
   return single_linkage_imatrix(file, seqnames, dconv, m);
 }

//' @export
 // [[Rcpp::export]]
 Rcpp::IntegerMatrix single_linkage_imatrix_cached(
     const std::string file,
     const Rcpp::CharacterVector &seqnames,
     const std::vector<double> &thresholds,
     const double precision
 ) {
   const CachedDistanceConverter dconv(thresholds, precision);
   const int m = thresholds.size();
   return single_linkage_matrix2(file, seqnames, dconv, m);
 }
#endif
