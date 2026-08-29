// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_CLUSTERALGORITHMFACTORY_H_INCLUDED
#define OPTIMOTU_CLUSTERALGORITHMFACTORY_H_INCLUDED

#include "ClusterAlgorithm.h"
#include <memory>
#include <cstddef>

enum class ClusterInstanceRole
{
  Leaf,
  Parent
};

class ClusterAlgorithmFactory {
protected:
  ClusterAlgorithmFactory(const DistanceConverter & dconv);
public:
  const DistanceConverter & dconv;

  virtual ~ClusterAlgorithmFactory() = default;

  virtual std::unique_ptr<SingleClusterAlgorithm> create(j_t n, int verbose = 0) const = 0;
  virtual std::unique_ptr<SingleClusterAlgorithm> create(init_matrix_t & im, int verbose = 0) const = 0;
  virtual std::size_t estimate_bytes(
      j_t n,
      ClusterInstanceRole role = ClusterInstanceRole::Parent) const;
};

class ClusterMatrixFactory : public ClusterAlgorithmFactory{
private:
  const bool binary_search;
  const int fill_type;
public:
  ClusterMatrixFactory(
    const DistanceConverter & dconv,
    const bool binary_search,
    const int fill_type
  );

  std::unique_ptr<SingleClusterAlgorithm> create(j_t n, int verbose = 0) const override;
  std::unique_ptr<SingleClusterAlgorithm> create(init_matrix_t & im, int verbose = 0) const override;
  std::size_t estimate_bytes(
      j_t n,
      ClusterInstanceRole role = ClusterInstanceRole::Parent) const override;
};

class ClusterIndexedMatrixFactory : public ClusterAlgorithmFactory{
public:
  ClusterIndexedMatrixFactory(const DistanceConverter & dconv);

  std::unique_ptr<SingleClusterAlgorithm> create(j_t n, int verbose = 0) const override;
  std::unique_ptr<SingleClusterAlgorithm> create(init_matrix_t & im, int verbose = 0) const override;
  std::size_t estimate_bytes(
      j_t n,
      ClusterInstanceRole role = ClusterInstanceRole::Parent) const override;
};

class ClusterTreeFactory : public ClusterAlgorithmFactory{
  const int test;
public:
  ClusterTreeFactory(const DistanceConverter & dconv, int test);

  std::unique_ptr<SingleClusterAlgorithm> create(j_t n, int verbose = 0) const override;
  std::unique_ptr<SingleClusterAlgorithm> create(init_matrix_t & im, int verbose = 0) const override;
  std::size_t estimate_bytes(
      j_t n,
      ClusterInstanceRole role = ClusterInstanceRole::Parent) const override;
};

class ClusterSLINKFactory : public ClusterAlgorithmFactory{
public:
  ClusterSLINKFactory(const DistanceConverter & dconv);

  std::unique_ptr<SingleClusterAlgorithm> create(j_t n, int verbose = 0) const override;
  std::unique_ptr<SingleClusterAlgorithm> create(init_matrix_t & im, int verbose = 0) const override;
  std::size_t estimate_bytes(
      j_t n,
      ClusterInstanceRole role = ClusterInstanceRole::Parent) const override;
};

#endif //OPTIMOTU_CLUSTERALGORITHMFACTORY_H_INCLUDED
