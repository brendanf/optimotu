#include "ClusterAlgorithmFactory.h"
#include "ClusterMatrix.h"
#include "ClusterIndexedMatrix.h"
#include "ClusterTree.h"

using CAF = ClusterAlgorithmFactory;
CAF::ClusterAlgorithmFactory(const DistanceConverter & dconv) :
  dconv{dconv} {}

using CMF = ClusterMatrixFactory;
CMF::ClusterMatrixFactory(
  const DistanceConverter & dconv,
  const bool binary_search,
  const int fill_type
) :
  CAF{dconv},
  binary_search{binary_search},
  fill_type{fill_type}
{}

std::unique_ptr<SingleClusterAlgorithm> CMF::create(j_t n) const {
  if (this->binary_search) {
    switch (this->fill_type) {
    case LINEAR_FILL:
      return std::make_unique<ClusterMatrix<true, LINEAR_FILL>>(dconv, n);
    case BINARY_FILL:
      return std::make_unique<ClusterMatrix<true, BINARY_FILL>>(this->dconv, n);
    case TOPDOWN_FILL:
      return std::make_unique<ClusterMatrix<true, TOPDOWN_FILL>>(this->dconv, n);
    default:
      Rcpp::stop("unknown fill type");
    }
  } else {
    switch (fill_type) {
    case LINEAR_FILL:
      return std::make_unique<ClusterMatrix<false, LINEAR_FILL>>(this->dconv, n);
    case BINARY_FILL:
      return std::make_unique<ClusterMatrix<false, BINARY_FILL>>(this->dconv, n);
    case TOPDOWN_FILL:
      return std::make_unique<ClusterMatrix<false, TOPDOWN_FILL>>(this->dconv, n);
    default:
      Rcpp::stop("unknown fill type");
    }
  }
}

std::unique_ptr<SingleClusterAlgorithm> CMF::create(init_matrix_t & im) const {
  if (this->binary_search) {
    switch (this->fill_type) {
    case LINEAR_FILL:
      return std::make_unique<ClusterMatrix<true, LINEAR_FILL, internal_matrix_t>>(this->dconv, im);
    case BINARY_FILL:
      return std::make_unique<ClusterMatrix<true, BINARY_FILL, internal_matrix_t>>(this->dconv, im);
    case TOPDOWN_FILL:
      return std::make_unique<ClusterMatrix<true, TOPDOWN_FILL, internal_matrix_t>>(this->dconv, im);
    default:
      Rcpp::stop("unknown fill type");
    }
  } else {
    switch (fill_type) {
    case LINEAR_FILL:
      return std::make_unique<ClusterMatrix<false, LINEAR_FILL, internal_matrix_t>>(this->dconv, im);
    case BINARY_FILL:
      return std::make_unique<ClusterMatrix<false, BINARY_FILL, internal_matrix_t>>(this->dconv, im);
    case TOPDOWN_FILL:
      return std::make_unique<ClusterMatrix<false, TOPDOWN_FILL, internal_matrix_t>>(this->dconv, im);
    default:
      Rcpp::stop("unknown fill type");
    }
  }
}

using CIMF = ClusterIndexedMatrixFactory;

CIMF::ClusterIndexedMatrixFactory(const DistanceConverter & dconv) : CAF{dconv} {}

std::unique_ptr<SingleClusterAlgorithm> CIMF::create(j_t n) const {
  return std::make_unique<ClusterIndexedMatrix<>>(this->dconv, n);
}

std::unique_ptr<SingleClusterAlgorithm> CIMF::create(init_matrix_t & im) const {
  return std::make_unique<ClusterIndexedMatrix<internal_matrix_t>>(this->dconv, im);
}

using CTF = ClusterTreeFactory;

CTF::ClusterTreeFactory(const DistanceConverter & dconv) : CAF{dconv} {}

std::unique_ptr<SingleClusterAlgorithm> CTF::create(j_t n) const {
  return std::make_unique<ClusterTree>(this->dconv, n);
}

std::unique_ptr<SingleClusterAlgorithm> CTF::create(init_matrix_t & im) const {
  return std::make_unique<ClusterTree>(this->dconv, im);
}


