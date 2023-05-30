#include "config.h"

double element_as_double(Rcpp::List obj, std::string e, std::string name) {
  if (obj[e] == R_NilValue) {
    Rcpp::stop(
      "invalid `" + name + "`: missing member `" + e + "`"
    );
  }
  if (!Rcpp::is<Rcpp::NumericVector>(obj[e])) {
    Rcpp::stop(
      "invalid `" + name + "`: member `" + e + "` is not numeric"
    );
  }
  auto item = Rcpp::as<Rcpp::NumericVector>(obj);
  if (item.length() != 1) {
    Rcpp::stop(
      "invalid `" + name + "`: length of member `" + e + "` is not 1"
    );
  }
  return item[0];
}

std::vector<double> element_as_double_vector(Rcpp::List obj, std::string e, std::string name) {
  if (obj[e] == R_NilValue) {
    Rcpp::stop(
      "invalid `" + name + "`: missing member `" + e + "`"
    );
  }
  if (!Rcpp::is<Rcpp::NumericVector>(obj[e])) {
    Rcpp::stop(
      "invalid `" + name + "`: member `" + e + "` is not numeric"
    );
  }
  return Rcpp::as<std::vector<double>>(obj);
}

int element_as_int(Rcpp::List obj, std::string e, std::string name) {
  if (obj[e] == R_NilValue) {
    Rcpp::stop(
      "invalid `" + name + "`: missing member `" + e + "`"
    );
  }
  if (!Rcpp::is<Rcpp::IntegerVector>(obj[e])) {
    Rcpp::stop(
      "invalid `" + name + "`: member `" + e + "` is not an integer"
    );
  }
  auto item = Rcpp::as<Rcpp::IntegerVector>(obj);
  if (item.length() != 1) {
    Rcpp::stop(
      "invalid `" + name + "`: length of member `" + e + "` is not 1"
    );
  }
  return item[0];
}

std::string element_as_string(Rcpp::List obj, std::string e, std::string name) {
  if (obj[e] == R_NilValue) {
    Rcpp::stop(
      "invalid `" + name + "`: missing member `" + e + "`"
    );
  }
  if (!Rcpp::is<Rcpp::CharacterVector>(obj[e])) {
    Rcpp::stop(
      "invalid `" + name + "`: member `" + e + "` is not of type character"
    );
  }
  auto item = Rcpp::as<Rcpp::CharacterVector>(obj);
  if (item.length() != 1) {
    Rcpp::stop(
      "invalid `" + name + "`: length of member `" + e + "` is not 1"
    );
  }
  return Rcpp::as<std::string>(item[0]);
}

bool element_as_bool(Rcpp::List obj, std::string e, std::string name) {
  if (obj[e] == R_NilValue) {
    Rcpp::stop(
      "invalid `" + name + "`: missing member `" + e + "`"
    );
  }
  if (!Rcpp::is<Rcpp::LogicalVector>(obj[e])) {
    Rcpp::stop(
      "invalid `" + name + "`: member `" + e + "` is not logical"
    );
  }
  auto item = Rcpp::as<Rcpp::LogicalVector>(obj);
  if (item.length() != 1) {
    Rcpp::stop(
      "invalid `" + name + "`: length of member `" + e + "` is not 1"
    );
  }
  return item[0];
}

std::unique_ptr<DistanceConverter> create_distance_converter(Rcpp::List config) {
  if (!config.inherits("optimotu_threshold_config")) {
    Rcpp::stop(
      "'threshold_config' must be of class 'optimotu_threshold_config'"
    );
  }
  std::string method = element_as_string(config, "method", "threshold_config");
  if (method == "uniform") {
    double from = element_as_double(config, "from", "threshold_uniform");
    double to = element_as_double(config, "to", "threshold_uniform");
    double by = element_as_double(config, "by", "threshold_uniform");
    return std::make_unique<UniformDistanceConverter>(from, to, by);
  }
  if (method == "set") {
    std::vector<double> thresholds =
      element_as_double_vector(config, "thresholds", "threshold_set");
    return std::make_unique<ArrayDistanceConverter>(thresholds);
  }
  if (method == "lookup") {
    std::vector<double> thresholds =
      element_as_double_vector(config, "thresholds", "threshold_lookup");
    double precision = element_as_double(config, "precision", "threshold_lookup");
    return std::make_unique<CachedDistanceConverter>(thresholds, precision);
  }
  Rcpp::stop("invalid `threshold_config`: unknown method: " + method);
}

std::unique_ptr<ClusterAlgorithmFactory> create_cluster_algorithm(
    Rcpp::List config,
    DistanceConverter * dconv
) {
  if (!config.inherits("optimotu_cluster_config")) {
    Rcpp::stop(
      "'cluster_config' must be of class 'optimotu_cluster_config'"
    );
  }
  std::string method = element_as_string(config, "method", "cluster_config");
  if (method == "matrix") {
    bool do_binary_search = element_as_bool(config, "binary_search", "cluster_matrix");
    int fill_type = element_as_int(config, "fill_method", "cluster_matrix");
    return std::make_unique<ClusterMatrixFactory>(*dconv, do_binary_search, fill_type);
  } else if (method == "index") {
    return std::make_unique<ClusterIndexedMatrixFactory>(*dconv);
  } else if (method == "tree") {
    return std::make_unique<ClusterTreeFactory>(*dconv);
  } else {
    Rcpp::stop("unknown cluster method");
  }
}

MultipleClusterAlgorithm create_multiple_cluster_algorithm(
  Rcpp::List parallel_config,
  ClusterAlgorithmFactory &factory,
  Rcpp::CharacterVector seqnames,
  Rcpp::ListOf<Rcpp::CharacterVector> subset_names
) {
  int threads = element_as_int(parallel_config, "threads", "parallel_config");
  return MultipleClusterAlgorithm(
    factory,
    Rcpp::as<std::vector<std::string>>(seqnames),
    Rcpp::as<std::vector<std::vector<std::string>>>(subset_names),
    threads
  );
}

std::unique_ptr<ClusterWorker> create_cluster_worker(
    Rcpp::List config,
    ClusterAlgorithm * algo,
    std::istream &file
) {
  if (!config.inherits("optimotu_parallel_config")) {
    Rcpp::stop(
      "'parallel_config' must be of class 'optimotu_parallel_config'"
    );
  }
  std::string method = element_as_string(config, "method", "parallel_config");
  int threads = element_as_int(config, "threads", "parallel_config");
  if (method == "merge") {
    return std::make_unique<MergeClusterWorker>(algo, file, threads);
  } else if (method == "concurrent") {
    return std::make_unique<ConcurrentClusterWorker>(algo, file, threads);
  } else if (method == "hierarchical") {
    int shards = element_as_int(config, "shards", "parallel_config");
    return std::make_unique<HierarchicalClusterWorker>(algo, file, threads, shards);
  } else {
    Rcpp::stop("unknown parallelization method");
  }
}
