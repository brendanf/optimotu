#ifndef OPTIMOTU_CONFIG_H_INCLUDED
#define OPTIMOTU_CONFIG_H_INCLUDED

#ifdef OPTIMOTU_R

#include <Rcpp.h>

#include "ClusterAlgorithmFactory.h"
#include "MultipleClusterAlgorithm.h"
#include "ClusterWorker.h"

#include "DistanceConverter.h"

std::unique_ptr<DistanceConverter> create_distance_converter(Rcpp::List config);
std::unique_ptr<ClusterAlgorithmFactory> create_cluster_algorithm(
    Rcpp::List config,
    DistanceConverter * dconv
);
std::unique_ptr<MultipleClusterAlgorithm> create_multiple_cluster_algorithm(
    Rcpp::List parallel_config,
    ClusterAlgorithmFactory &factory,
    Rcpp::CharacterVector seqnames,
    Rcpp::ListOf<Rcpp::CharacterVector> subset_names
);
std::unique_ptr<ClusterWorker> create_cluster_worker(
  Rcpp::List config,
  ClusterAlgorithm * algo,
  std::istream &file
);

#endif //OPTIMOTU_R

#endif //OPTIMOTU_CONFIG_H_INCLUDED
