#ifndef OPTIMOTU_CONFIG_H_INCLUDED
#define OPTIMOTU_CONFIG_H_INCLUDED

#include "ClusterAlgorithmFactory.h"
#include "MultipleClusterAlgorithm.h"
#include "ClusterWorker.h"

#include "DistanceConverter.h"

#include "AlignClusterWorker.h"

#ifdef OPTIMOTU_R
#include <Rcpp.h>
#endif //OPTIMOTU_R

std::unique_ptr<DistanceConverter> create_distance_converter(
    const std::string &type,
    const double from,
    const double to,
    const double by,
    std::vector<double> thresholds,
    double precision
);
std::unique_ptr<ClusterAlgorithmFactory> create_cluster_algorithm(
    const std::string &method,
    const bool do_binary_search,
    const int fill_type,
    DistanceConverter * dconv
);
std::unique_ptr<MultipleClusterAlgorithm> create_multiple_cluster_algorithm(
    const int threads,
    ClusterAlgorithmFactory &factory,
    const std::vector<std::string> seqnames,
    const std::vector<std::vector<std::string>> subset_names
);
std::unique_ptr<ClusterWorker> create_cluster_worker(
    const std::string &method,
    const int threads,
    const int shards,
    ClusterAlgorithm * algo,
    std::istream &file
);
std::unique_ptr<AlignClusterWorker> create_align_cluster_worker(
    const std::string &type,
    const std::vector<std::string> &seq,
    const double breakpoint,
    SingleClusterAlgorithm &cluster,
    const uint8_t threads
);

#ifdef OPTIMOTU_R
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
std::unique_ptr<AlignClusterWorker> create_align_cluster_worker(
    Rcpp::List dist_config,
    Rcpp::List parallel_config,
    const std::vector<std::string> &seq,
    SingleClusterAlgorithm &cluster,
    bool verbose = false
);

#endif //OPTIMOTU_R

#endif //OPTIMOTU_CONFIG_H_INCLUDED
