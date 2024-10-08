#ifdef OPTIMOTU_R

#include <fstream>
#include <Rcpp.h>
#include "single_linkage.h"
#include "config.h"

// [[Rcpp::export]]
Rcpp::RObject distmx_cluster_single(
  const std::string file,
  const Rcpp::CharacterVector seqnames,
  const Rcpp::List threshold_config,
  const Rcpp::List clust_config,
  const Rcpp::List parallel_config,
  const std::string output_type = "matrix",
  const bool verbose = false
) {
  if (output_type != "matrix" && output_type != "hclust") {
    OPTIMOTU_STOP("Unknown 'output_type'");
  }
  std::ifstream infile(file);
  if (!infile.is_open()) {
    OPTIMOTU_STOP("failed to open input file");
  }
  auto dconv = create_distance_converter(threshold_config);
  Rcpp::IntegerMatrix im(dconv->m, seqnames.size());
  auto algo = create_cluster_algorithm(clust_config, dconv.get())->create(im);
  auto worker = create_cluster_worker(parallel_config, algo.get(), infile);
  int threads = worker->n_threads();
  if (threads == 1) {
    (*worker)(0, 1);
  } else {
    RcppParallel::parallelFor(0, threads, *worker, 1, threads);
  }
  worker->finalize();
  Rcpp::RObject output = R_NilValue;
  if (output_type == "matrix") {
    auto method = (Rcpp::as<Rcpp::CharacterVector>(clust_config["method"]))[0];
    if (method == "tree" || method == "slink") {
      internal_matrix_t m(im);
      algo->write_to_matrix(m);
    }
    output = im;
  } else if (output_type == "hclust") {
    output = algo->as_hclust(seqnames);
  }
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
    const bool verbose = false
) {
  if (output_type != "matrix" && output_type != "hclust") {
    OPTIMOTU_STOP("Unknown 'output_type'");
  }
  Rcpp::List output;
  std::ifstream infile(file);
  if (!infile.is_open()) {
    OPTIMOTU_STOP("failed to open input file");
  }
  // OPTIMOTU_COUT << "creating DistanceConverter..." << std::flush;
  auto dconv = create_distance_converter(threshold_config);
  // OPTIMOTU_COUT << "done" << std::endl
              // << "creating ClusterAlgorithmFactory..." << std::flush;
  auto factory = create_cluster_algorithm(method_config, dconv.get());
  // OPTIMOTU_COUT << "done" << std::endl
  //             << "creating MultipleClusterAlgorithm..." << std::flush;
  auto algo = create_multiple_cluster_algorithm(parallel_config, *factory, seqnames, which);
  // OPTIMOTU_COUT << "done" << std::endl
  //             << "creating ClusterWorker..." << std::flush;
  auto worker = create_cluster_worker(parallel_config, algo.get(), infile);
  // OPTIMOTU_COUT << "done" << std::endl
  //             << "clustering..." << std::flush;

  if (worker->n_threads() == 1) {
    (*worker)(0, 1);
  } else {
    RcppParallel::parallelFor(0, worker->n_threads(), *worker, 1, worker->n_threads());
  }
  // OPTIMOTU_COUT << "done" << std::endl
  //             << "finalizing worker..." << std::flush;
  worker->finalize();
  // OPTIMOTU_COUT << "done" << std::endl
  //             << "creating output..." << std::flush;
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
