#ifdef OPTIMOTU_R
#include <fstream>
#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include "ClusterMatrix.h"
#include "ClusterIndexedMatrix.h"

typedef RcppParallel::RMatrix<int> matrix_t;

std::unique_ptr<ClusterAlgorithm> create_clustermatrix(
    const DistanceConverter &dconv,
    Rcpp::IntegerMatrix im,
    const bool do_binary_search = true,
    const int fill_type = LINEAR_FILL
) {
  if (do_binary_search) {
    switch (fill_type) {
    case LINEAR_FILL:
      return std::make_unique<ClusterMatrix<matrix_t, true, LINEAR_FILL>>(dconv, im);
    case BINARY_FILL:
      return std::make_unique<ClusterMatrix<matrix_t, true, BINARY_FILL>>(dconv, im);
    case TOPDOWN_FILL:
      return std::make_unique<ClusterMatrix<matrix_t, true, TOPDOWN_FILL>>(dconv, im);
    default:
      Rcpp::stop("unknown fill type");
    }
  } else {
    switch (fill_type) {
    case LINEAR_FILL:
      return std::make_unique<ClusterMatrix<matrix_t, false, LINEAR_FILL>>(dconv, im);
    case BINARY_FILL:
      return std::make_unique<ClusterMatrix<matrix_t, false, BINARY_FILL>>(dconv, im);
    case TOPDOWN_FILL:
      return std::make_unique<ClusterMatrix<matrix_t, false, TOPDOWN_FILL>>(dconv, im);
    default:
      Rcpp::stop("unknown fill type");
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

Rcpp::IntegerMatrix distmx_cluster_matrix2(
    const std::string file,
    const Rcpp::CharacterVector &seqnames,
    const DistanceConverter &dconv,
    const int m,
    const int threads = 1,
    const bool do_binary_search = true,
    const int fill_method = BINARY_FILL
) {
  Rcpp::IntegerMatrix im(m, seqnames.size());
  auto cm = create_clustermatrix(dconv, im, do_binary_search, fill_method);
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
//' @rdname distmx_cluster
// [[Rcpp::export]]
Rcpp::IntegerMatrix distmx_cluster_matrix2_uniform(
      const std::string file,
      const Rcpp::CharacterVector &seqnames,
      const float dmin,
      const float dmax,
      const float dstep,
      const int threads = 1,
      const bool do_binary_search = true,
      const int fill_method = 2
) {
   const UniformDistanceConverter dconv(dmin, dmax, dstep);
   const int m = (int) ceil((dmax - dmin)/dstep) + 1;
   return distmx_cluster_matrix2(file, seqnames, dconv, m, threads,
                                 do_binary_search, fill_method);
}

//' @export
//' @rdname distmx_cluster
// [[Rcpp::export]]
Rcpp::IntegerMatrix distmx_cluster_matrix2_array(
      const std::string file,
      const Rcpp::CharacterVector &seqnames,
      const std::vector<double> &thresholds,
      const int threads = 1,
      const bool do_binary_search = true,
      const int fill_method = 2
) {
   const ArrayDistanceConverter dconv(thresholds);
   const int m = thresholds.size();
   return distmx_cluster_matrix2(file, seqnames, dconv, m, threads,
                                 do_binary_search, fill_method);
}

//' @export
//' @rdname distmx_cluster
// [[Rcpp::export]]
Rcpp::IntegerMatrix distmx_cluster_matrix2_cached(
      const std::string file,
      const Rcpp::CharacterVector &seqnames,
      const std::vector<double> &thresholds,
      const double precision,
      const int threads = 1,
      const bool do_binary_search = true,
      const int fill_method = 2
) {
   const CachedDistanceConverter dconv(thresholds, precision);
   const int m = thresholds.size();
   return distmx_cluster_matrix2(file, seqnames, dconv, m, threads,
                                 do_binary_search, fill_method);
}

Rcpp::IntegerMatrix distmx_cluster_imatrix(
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
//' @rdname distmx_cluster
// [[Rcpp::export]]
Rcpp::IntegerMatrix distmx_cluster_imatrix_uniform(
    const std::string file,
    const Rcpp::CharacterVector &seqnames,
    const float dmin,
    const float dmax,
    const float dstep
) {
  const UniformDistanceConverter dconv(dmin, dmax, dstep);
  const int m = (int) ceilf((dmax - dmin)/dstep) + 1;
  return distmx_cluster_imatrix(file, seqnames, dconv, m);
}

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix distmx_cluster_imatrix_array(
    const std::string file,
    const Rcpp::CharacterVector &seqnames,
    const std::vector<double> &thresholds
) {
  const ArrayDistanceConverter dconv(thresholds);
  const int m = thresholds.size();
  return distmx_cluster_imatrix(file, seqnames, dconv, m);
}

//' @export
//' @rdname distmx_cluster
// [[Rcpp::export]]
Rcpp::IntegerMatrix distmx_cluster_imatrix_cached(
    const std::string file,
    const Rcpp::CharacterVector &seqnames,
    const std::vector<double> &thresholds,
    const double precision
) {
  const CachedDistanceConverter dconv(thresholds, precision);
  const int m = thresholds.size();
  return distmx_cluster_matrix2(file, seqnames, dconv, m);
}
#endif
