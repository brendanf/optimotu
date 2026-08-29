// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "ClusterAlgorithmFactory.h"
#include "ClusterMatrix.h"
#include "ClusterIndexedMatrix.h"
#include "ClusterTree.h"
#include "ClusterSLINK.h"
#include <cstdint>

using CAF = ClusterAlgorithmFactory;
CAF::ClusterAlgorithmFactory(const DistanceConverter & dconv) :
  dconv{dconv} {}

std::size_t CAF::estimate_bytes(j_t n, ClusterInstanceRole) const
{
  const std::size_t ns = static_cast<std::size_t>(n);
  const std::size_t ms = static_cast<std::size_t>(dconv.m);
  return ns * ms * sizeof(int);
}

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

std::size_t CMF::estimate_bytes(j_t n, ClusterInstanceRole) const
{
  const std::size_t ns = static_cast<std::size_t>(n);
  const std::size_t ms = static_cast<std::size_t>(dconv.m);
  // Matrix plus helper vectors and some allocator overhead.
  return ns * ms * sizeof(int) + ns * sizeof(int) * 8;
}

using CIMF = ClusterIndexedMatrixFactory;

CIMF::ClusterIndexedMatrixFactory(const DistanceConverter & dconv) : CAF{dconv} {}

std::unique_ptr<SingleClusterAlgorithm> CIMF::create(j_t n, int verbose) const {
  return create_cluster_indexed_matrix(this->dconv, n, verbose);
}

std::unique_ptr<SingleClusterAlgorithm> CIMF::create(init_matrix_t & im, int) const {
  return create_cluster_indexed_matrix(this->dconv, im);
}

std::size_t CIMF::estimate_bytes(j_t n, ClusterInstanceRole) const
{
  const std::size_t ns = static_cast<std::size_t>(n);
  const std::size_t ms = static_cast<std::size_t>(dconv.m);
  // Indexed matrix stores matrix + index list + buffers.
  return ns * ms * sizeof(int) + ns * sizeof(int) * 12;
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

std::size_t CTF::estimate_bytes(j_t n, ClusterInstanceRole) const
{
  const std::size_t ns = static_cast<std::size_t>(n);
  // ClusterTree uses roughly 2*n nodes and free-list structures.
  return ns * sizeof(std::uint64_t) * 32;
}

using CSF = ClusterSLINKFactory;

CSF::ClusterSLINKFactory(const DistanceConverter & dconv) : CAF{dconv} {}

std::unique_ptr<SingleClusterAlgorithm> CSF::create(j_t n, int verbose) const {
  return create_cluster_slink(this->dconv, n, verbose);
}

std::unique_ptr<SingleClusterAlgorithm> CSF::create(init_matrix_t & im, int verbose) const {
  return create_cluster_slink(this->dconv, im, verbose);
}

std::size_t CSF::estimate_bytes(j_t n, ClusterInstanceRole role) const
{
  const std::size_t ns = static_cast<std::size_t>(n);
  if (role == ClusterInstanceRole::Leaf)
  {
    // Pi, Lambda, M vectors plus allocator overhead.
    return ns * sizeof(std::uint64_t) * 4;
  }
  // Delegate tree storage for parent instances.
  return ns * sizeof(std::uint64_t) * 32;
}
