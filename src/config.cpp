#include "config.h"
#include "Wfa2ClusterWorker.h"
#include "EdlibClusterWorker.h"
#include "HybridClusterWorker.h"
#include "HammingClusterWorker.h"

std::unique_ptr<DistanceConverter> create_distance_converter(
    const std::string &type,
    const double from,
    const double to,
    const double by,
    std::vector<double> thresholds,
    double precision
) {
  if (type == "uniform") {
    return std::make_unique<UniformDistanceConverter>(from, to, by);
  }
  if (type == "set") {
    return std::make_unique<ArrayDistanceConverter>(thresholds);
  }
  if (type == "lookup") {
    return std::make_unique<CachedDistanceConverter>(thresholds, precision);
  }
  OPTIMOTU_STOP("invalid `threshold_config`: unknown type: " + type);
}

std::unique_ptr<ClusterAlgorithmFactory> create_cluster_algorithm(
    const std::string &method,
    const bool do_binary_search,
    const int fill_type,
    DistanceConverter * dconv
) {
  if (method == "matrix") {
    return std::make_unique<ClusterMatrixFactory>(*dconv, do_binary_search, fill_type);
  } else if (method == "index") {
    return std::make_unique<ClusterIndexedMatrixFactory>(*dconv);
  } else if (method == "tree") {
    return std::make_unique<ClusterTreeFactory>(*dconv);
  } else if (method == "slink") {
    return std::make_unique<ClusterSLINKFactory>(*dconv);
  } else {
    OPTIMOTU_STOP("unknown cluster method");
  }
}

std::unique_ptr<MultipleClusterAlgorithm> create_multiple_cluster_algorithm(
    const int threads,
    ClusterAlgorithmFactory &factory,
    const std::vector<std::string> seqnames,
    const std::vector<std::vector<std::string>> subset_names
) {
  auto out = std::make_unique<MultipleClusterAlgorithm>(
    factory,
    seqnames,
    subset_names,
    threads
  );
  return out;
}

std::unique_ptr<ClusterWorker> create_cluster_worker(
    const std::string &method,
    const int threads,
    const int shards,
    ClusterAlgorithm * algo,
    std::istream &file
) {
  if (method == "merge") {
    return std::make_unique<MergeClusterWorker>(algo, file, threads);
  } else if (method == "concurrent") {
    return std::make_unique<ConcurrentClusterWorker>(algo, file, threads);
  } else if (method == "hierarchical") {
    return std::make_unique<HierarchicalClusterWorker>(algo, file, threads, shards);
  } else {
    OPTIMOTU_STOP("unknown parallelization method");
  }
}

std::unique_ptr<AlignClusterWorker> create_align_cluster_worker(
    const std::string &type,
    const std::vector<std::string> &seq,
    const double breakpoint,
    ClusterAlgorithm &cluster,
    const std::uint8_t threads
) {
  if (type == "split") {
    return std::make_unique<HybridSplitClusterWorker>(seq, cluster, threads, breakpoint);
  } else if (type == "concurrent") {
    return std::make_unique<HybridConcurrentClusterWorker>(seq, cluster, threads, breakpoint);
  } else {
    Rcpp::stop("invalid parallel type");
  }
}

#ifdef OPTIMOTU_R
double element_as_double(Rcpp::List obj, std::string e, std::string name) {
  if (!obj.containsElementNamed(e.c_str())) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: missing member `" + e + "`"
    );
  }
  Rcpp::RObject item = obj[e];
  if (!Rcpp::is<Rcpp::NumericVector>(item)) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: member `" + e + "` is not numeric"
    );
  }
  auto item_dbl = Rcpp::as<Rcpp::NumericVector>(item);
  if (item_dbl.length() != 1) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: length of member `" + e + "` is not 1"
    );
  }
  return item_dbl[0];
}

std::vector<double> element_as_double_vector(Rcpp::List obj, std::string e, std::string name) {
  if (!obj.containsElementNamed(e.c_str())) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: missing member `" + e + "`"
    );
  }
  Rcpp::RObject item = obj[e];
  if (!Rcpp::is<Rcpp::NumericVector>(item)) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: member `" + e + "` is not numeric"
    );
  }
  return Rcpp::as<std::vector<double>>(item);
}

int element_as_int(Rcpp::List obj, std::string e, std::string name) {
  if (!obj.containsElementNamed(e.c_str())) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: missing member `" + e + "`"
    );
  }
  Rcpp::RObject item = obj[e];
  if (!Rcpp::is<Rcpp::IntegerVector>(item)) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: member `" + e + "` is not an integer"
    );
  }
  auto item_int = Rcpp::as<Rcpp::IntegerVector>(item);
  if (item_int.length() != 1) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: length of member `" + e + "` is not 1"
    );
  }
  return item_int[0];
}

std::string element_as_string(Rcpp::List obj, std::string e, std::string name) {
  if (!obj.containsElementNamed(e.c_str())) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: missing member `" + e + "`"
    );
  }
  Rcpp::RObject item = obj[e];
  if (!Rcpp::is<Rcpp::CharacterVector>(item)) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: member `" + e + "` is not of type character"
    );
  }
  auto item_str = Rcpp::as<Rcpp::CharacterVector>(item);
  if (item_str.length() != 1) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: length of member `" + e + "` is not 1"
    );
  }
  return Rcpp::as<std::string>(item_str[0]);
}

bool element_as_bool(Rcpp::List obj, std::string e, std::string name) {
  if (!obj.containsElementNamed(e.c_str())) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: missing member `" + e + "`"
    );
  }
  Rcpp::RObject item = obj[e];
  if (!Rcpp::is<Rcpp::LogicalVector>(item)) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: member `" + e + "` is not logical"
    );
  }
  auto item_bool = Rcpp::as<Rcpp::LogicalVector>(item);
  if (item_bool.length() != 1) {
    OPTIMOTU_STOP(
      "invalid `" + name + "`: length of member `" + e + "` is not 1"
    );
  }
  return item_bool[0];
}

std::unique_ptr<DistanceConverter> create_distance_converter(Rcpp::List config) {
  if (!config.inherits("optimotu_threshold_config")) {
    OPTIMOTU_STOP(
      "'threshold_config' must be of class 'optimotu_threshold_config'"
    );
  }
  std::string type = element_as_string(config, "type", "threshold_config");
  if (type == "uniform") {
    double from = element_as_double(config, "from", "threshold_uniform");
    double to = element_as_double(config, "to", "threshold_uniform");
    double by = element_as_double(config, "by", "threshold_uniform");
    return std::make_unique<UniformDistanceConverter>(from, to, by);
  }
  if (type == "set") {
    std::vector<double> thresholds =
      element_as_double_vector(config, "thresholds", "threshold_set");
    return std::make_unique<ArrayDistanceConverter>(thresholds);
  }
  if (type == "lookup") {
    std::vector<double> thresholds =
      element_as_double_vector(config, "thresholds", "threshold_lookup");
    double precision = element_as_double(config, "precision", "threshold_lookup");
    return std::make_unique<CachedDistanceConverter>(thresholds, precision);
  }
  OPTIMOTU_STOP("invalid `threshold_config`: unknown type: " + type);
}

std::unique_ptr<ClusterAlgorithmFactory> create_cluster_algorithm(
    Rcpp::List config,
    DistanceConverter * dconv
) {
  if (!config.inherits("optimotu_cluster_config")) {
    OPTIMOTU_STOP(
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
  } else if (method == "slink") {
    return std::make_unique<ClusterSLINKFactory>(*dconv);
  } else {
    OPTIMOTU_STOP("unknown cluster method");
  }
}

std::unique_ptr<MultipleClusterAlgorithm> create_multiple_cluster_algorithm(
  Rcpp::List parallel_config,
  ClusterAlgorithmFactory &factory,
  Rcpp::CharacterVector seqnames,
  Rcpp::ListOf<Rcpp::CharacterVector> subset_names
) {
  int threads = element_as_int(parallel_config, "threads", "parallel_config");
  auto out = std::make_unique<MultipleClusterAlgorithm>(
    factory,
    Rcpp::as<std::vector<std::string>>(seqnames),
    Rcpp::as<std::vector<std::vector<std::string>>>(subset_names),
    threads
  );
  return out;
}

std::unique_ptr<ClusterWorker> create_cluster_worker(
    Rcpp::List config,
    ClusterAlgorithm * algo,
    std::istream &file
) {
  if (!config.inherits("optimotu_parallel_config")) {
    OPTIMOTU_STOP(
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
    OPTIMOTU_STOP("unknown parallelization method");
  }
}

std::unique_ptr<AlignClusterWorker> create_align_cluster_worker(
    Rcpp::List dist_config,
    Rcpp::List parallel_config,
    const std::vector<std::string> &seq,
    ClusterAlgorithm &cluster,
    bool verbose
) {
  if (!dist_config.inherits("optimotu_dist_config")) {
    OPTIMOTU_STOP(
      "'dist_config' must be of class 'optimotu_dist_config'"
    );
  }
  if (!parallel_config.inherits("optimotu_parallel_config")) {
    OPTIMOTU_STOP(
      "'parallel_config' must be of class 'optimotu_parallel_config'"
    );
  }
  std::string dist_method = element_as_string(dist_config, "method", "dist_config");
  std::string par_method = element_as_string(parallel_config, "method", "parallel_config");
  int threads = element_as_int(parallel_config, "threads", "parallel_config");
  if (dist_method == "wfa2") {
    int match = element_as_int(dist_config, "match", "dist_config");
    int mismatch = element_as_int(dist_config, "mismatch", "dist_config");
    int gap_open = element_as_int(dist_config, "gap_open", "dist_config");
    int gap_extend = element_as_int(dist_config, "gap_extend", "dist_config");
    int gap_open2 = element_as_int(dist_config, "gap_open2", "dist_config");
    int gap_extend2 = element_as_int(dist_config, "gap_extend2", "dist_config");
    if (par_method == "merge") {
      return std::make_unique<Wfa2SplitClusterWorker>(
        seq,
        cluster, threads, match, mismatch,
        gap_open, gap_extend, gap_open2, gap_extend2, verbose
      );
    } else if (par_method == "concurrent") {
      return std::make_unique<Wfa2ConcurrentClusterWorker>(
        seq,
        cluster, threads, match, mismatch,
        gap_open, gap_extend, gap_open2, gap_extend2, verbose
      );
    } else {
      OPTIMOTU_STOP("unknown parallelization method");
    }
  } else if (dist_method == "edlib") {
    if (par_method == "merge") {
      return std::make_unique<EdlibSplitClusterWorker>(
        seq,
        cluster, threads, verbose
      );
    } else if (par_method == "concurrent") {
      return std::make_unique<EdlibConcurrentClusterWorker>(
        seq,
        cluster, threads, verbose
      );
    } else {
      OPTIMOTU_STOP("unknown parallelization method");
    }
  } else if (dist_method == "hybrid") {
    double breakpoint = element_as_double(dist_config, "cutoff", "dist_config");
    if (par_method == "merge") {
      return std::make_unique<HybridSplitClusterWorker>(
        seq,
        cluster, threads, breakpoint, verbose
      );
    } else if (par_method == "concurrent") {
      return std::make_unique<HybridConcurrentClusterWorker>(
        seq,
        cluster, threads, breakpoint, verbose
      );
    } else {
      OPTIMOTU_STOP("unknown parallelization method");
    }
  } else if (dist_method == "hamming") {
    int min_overlap = element_as_int(dist_config, "min_overlap", "dist_config");
    bool ignore_gaps = element_as_bool(dist_config, "ignore_gaps", "dist_config");
    if (par_method == "merge") {
      return std::make_unique<HammingSplitClusterWorker>(
        seq,
        cluster, threads, min_overlap, ignore_gaps, verbose
      );
    } else if (par_method == "concurrent") {
      return std::make_unique<HammingConcurrentClusterWorker>(
        seq,
        cluster, threads, min_overlap, ignore_gaps, verbose
      );
    } else {
      OPTIMOTU_STOP("unknown parallelization method");
    }
  } else {
    OPTIMOTU_STOP("unknown sequence-cluster method");
  }
}
#endif // OPTIMOTU_R
