// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_CLUSTERALGORITHMFACTORY_H_INCLUDED
#define OPTIMOTU_CLUSTERALGORITHMFACTORY_H_INCLUDED

#include "ClusterAlgorithm.h"
#include <memory>

class ClusterAlgorithmFactory {
protected:
  ClusterAlgorithmFactory(const DistanceConverter & dconv);
public:
  const DistanceConverter & dconv;

  virtual ~ClusterAlgorithmFactory() = default;

  virtual std::unique_ptr<SingleClusterAlgorithm> create(j_t n) const = 0;
  virtual std::unique_ptr<SingleClusterAlgorithm> create(init_matrix_t & im) const = 0;
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

  std::unique_ptr<SingleClusterAlgorithm> create(j_t n) const override;
  std::unique_ptr<SingleClusterAlgorithm> create(init_matrix_t & im) const override;
};

class ClusterIndexedMatrixFactory : public ClusterAlgorithmFactory{
public:
  ClusterIndexedMatrixFactory(const DistanceConverter & dconv);

  std::unique_ptr<SingleClusterAlgorithm> create(j_t n) const override;
  std::unique_ptr<SingleClusterAlgorithm> create(init_matrix_t & im) const override;
};

class ClusterTreeFactory : public ClusterAlgorithmFactory{
  const int verbose;
  const int test;
public:
  ClusterTreeFactory(const DistanceConverter & dconv, int verbose, int test);

  std::unique_ptr<SingleClusterAlgorithm> create(j_t n) const override;
  std::unique_ptr<SingleClusterAlgorithm> create(init_matrix_t & im) const override;
};

class ClusterSLINKFactory : public ClusterAlgorithmFactory{
public:
  ClusterSLINKFactory(const DistanceConverter & dconv);

  std::unique_ptr<SingleClusterAlgorithm> create(j_t n) const override;
  std::unique_ptr<SingleClusterAlgorithm> create(init_matrix_t & im) const override;
};

#endif //OPTIMOTU_CLUSTERALGORITHMFACTORY_H_INCLUDED
