// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "ClusterAlgorithmFactory.h"
#include "ClusterMatrix.h"
#include "ClusterIndexedMatrix.h"
#include "ClusterTree.h"
#include "ClusterSLINK.h"

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

std::unique_ptr<SingleClusterAlgorithm> CMF::create(j_t n, int verbose) const {
  return create_cluster_matrix(this->dconv, n, this->binary_search, this->fill_type, verbose);
}

std::unique_ptr<SingleClusterAlgorithm> CMF::create(init_matrix_t & im, int) const {
  return create_cluster_matrix(this->dconv, im, this->binary_search, this->fill_type);
}

using CIMF = ClusterIndexedMatrixFactory;

CIMF::ClusterIndexedMatrixFactory(const DistanceConverter & dconv) : CAF{dconv} {}

std::unique_ptr<SingleClusterAlgorithm> CIMF::create(j_t n, int verbose) const {
  return create_cluster_indexed_matrix(this->dconv, n, verbose);
}

std::unique_ptr<SingleClusterAlgorithm> CIMF::create(init_matrix_t & im, int) const {
  return create_cluster_indexed_matrix(this->dconv, im);
}

using CTF = ClusterTreeFactory;

CTF::ClusterTreeFactory(const DistanceConverter & dconv, int test) :
  CAF{dconv}, test(test) {}

std::unique_ptr<SingleClusterAlgorithm> CTF::create(j_t n, int verbose) const {
  return create_cluster_tree(this->dconv, n, verbose, test);
}

std::unique_ptr<SingleClusterAlgorithm> CTF::create(init_matrix_t & im, int verbose) const {
  return create_cluster_tree(this->dconv, im, verbose, test);
}

using CSF = ClusterSLINKFactory;

CSF::ClusterSLINKFactory(const DistanceConverter & dconv) : CAF{dconv} {}

std::unique_ptr<SingleClusterAlgorithm> CSF::create(j_t n, int verbose) const {
  return create_cluster_slink(this->dconv, n, verbose);
}

std::unique_ptr<SingleClusterAlgorithm> CSF::create(init_matrix_t & im, int verbose) const {
  return create_cluster_slink(this->dconv, im, verbose);
}

