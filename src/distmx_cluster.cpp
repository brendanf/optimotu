#ifdef OPTIMOTU_R

#include <fstream>
#include <Rcpp.h>
#include "single_linkage.h"
#include "config.h"

// [[Rcpp::export]]
Rcpp::RObject distmx_cluster_single(
  const std::string file,
  Rcpp::CharacterVector seqnames,
  Rcpp::List threshold_config,
  Rcpp::List clust_config,
  Rcpp::List parallel_config,
  const std::string output_type = "matrix",
  const bool verbose = false,
  const bool by_name = false
) {
  if (output_type != "matrix" && output_type != "hclust") {
    OPTIMOTU_STOP("Unknown 'output_type'");
  }
  if (!clust_config.inherits("optimotu_cluster_config")) {
    OPTIMOTU_STOP("clust_config is not a valid cluster configuration");
  }
  if (!threshold_config.inherits("optimotu_threshold_config")) {
    OPTIMOTU_STOP("threshold_config is not a valid threshold configuration");
  }
  if (!parallel_config.inherits("optimotu_parallel_config")) {
    OPTIMOTU_STOP("parallel_config is not a valid parallel configuration");
  }
  if (!clust_config.containsElementNamed("method")) {
    OPTIMOTU_STOP("clust_config does not contain 'method' element");
  }
  if (!Rcpp::is<Rcpp::CharacterVector>(clust_config["method"])) {
    OPTIMOTU_STOP("clust_config$method is not a character vector");
  }
  std::ifstream infile(file);
  if (!infile.is_open()) {
    OPTIMOTU_STOP("failed to open input file");
  }
  if (verbose) OPTIMOTU_CERR << "creating DistanceConverter..." << std::flush;
  auto dconv = create_distance_converter(threshold_config);
  if (verbose) OPTIMOTU_CERR << "done\ncreating ClusterAlgorithm..." << std::flush;
  Rcpp::IntegerMatrix im(dconv->m, seqnames.size());
  auto algo = create_cluster_algorithm(clust_config, dconv.get())->create(im);
  if (verbose) OPTIMOTU_CERR << "done\ncreating ClusterWorker..." << std::flush;
  auto worker = create_cluster_worker<std::istream>(parallel_config, algo.get(), infile, by_name, seqnames);
  if (verbose) OPTIMOTU_CERR << "done\nclustering..." << std::flush;
  int threads = worker->n_threads();
  if (threads == 1) {
    (*worker)(0, 1);
  } else {
    RcppParallel::parallelFor(0, threads, *worker, 1, threads);
  }
  worker->finalize();
  if (verbose) OPTIMOTU_CERR << "done\ncreating output..." << std::flush;
  Rcpp::RObject output = R_NilValue;
  if (output_type == "matrix") {
    auto method = element_as_string(clust_config, "method", "clust_config");
    if (method == "tree" || method == "slink") {
      internal_matrix_t m(im);
      algo->write_to_matrix(m);
    }
    output = im;
  } else if (output_type == "hclust") {
    output = algo->as_hclust(seqnames);
  }
  if (verbose) OPTIMOTU_CERR << "done" << std::endl;
  return output;
}

// [[Rcpp::export]]
Rcpp::List distmx_cluster_multi(
    const std::string file,
    const Rcpp::CharacterVector seqnames,
    const Rcpp::ListOf<Rcpp::CharacterVector> which,
    const Rcpp::List threshold_config,
    const Rcpp::List method_config,
    const Rcpp::List parallel_config,
    const std::string output_type = "matrix",
    const bool verbose = false,
    const bool by_name = false
) {
  if (output_type != "matrix" && output_type != "hclust") {
    OPTIMOTU_STOP("Unknown 'output_type'");
  }
  Rcpp::List output;
  std::ifstream infile(file);
  if (!infile.is_open()) {
    OPTIMOTU_STOP("failed to open input file");
  }
  if (verbose) OPTIMOTU_CERR << "creating DistanceConverter..." << std::flush;
  auto dconv = create_distance_converter(threshold_config);
  if (verbose) OPTIMOTU_CERR << "done\ncreating ClusterAlgorithmFactory..."
                             << std::flush;
  auto factory = create_cluster_algorithm(method_config, dconv.get());
  if (verbose) OPTIMOTU_CERR << "done\ncreating MultipleClusterAlgorithm..."
                             << std::flush;
  auto algo = create_multiple_cluster_algorithm(parallel_config, *factory, seqnames, which);
  if (verbose) OPTIMOTU_CERR << "done\ncreating ClusterWorker..." << std::flush;
  auto worker = create_cluster_worker<std::istream>(parallel_config, algo.get(), infile, by_name, seqnames);
  if (verbose) OPTIMOTU_CERR << "done\nclustering..." << std::flush;
  if (worker->n_threads() == 1) {
    (*worker)(0, 1);
  } else {
    RcppParallel::parallelFor(0, worker->n_threads(), *worker, 1, worker->n_threads());
  }
  if (verbose) OPTIMOTU_CERR << "done\nfinalizing worker..." << std::flush;
  worker->finalize();
  if (verbose) OPTIMOTU_CERR << "done\ncreating output..." << std::flush;
  if (output_type == "matrix") {
    auto outlist = Rcpp::List(which.size());
    auto internal_out = std::vector<RcppParallel::RMatrix<int>>();
    internal_out.reserve(which.size());
    for (int i = 0; i < which.size(); ++i) {
      auto outi = Rcpp::IntegerMatrix(dconv->m, which[i].size());
      outlist[i] = outi;
      internal_out.emplace_back(outi);
    }
    algo->write_to_matrix(internal_out);
    output = outlist;
  } else if (output_type == "hclust") {
    output = algo->as_hclust();
  }
  output.names() = which.names();
  return output;
}
#endif // OPTIMOTU_R
