#include <fstream>
#include <Rcpp.h>
#include "single_linkage.h"
#include "config.h"

// [[Rcpp::export]]
Rcpp::RObject distmx_cluster_single(
  const std::string file,
  const Rcpp::CharacterVector seqnames,
  const Rcpp::List threshold_config,
  const Rcpp::List method_config,
  const Rcpp::List parallel_config,
  const std::string output_type = "matrix"
) {
  if (output_type != "matrix" && output_type != "hclust") {
    Rcpp::stop("Unknown 'output_type'");
  }
  std::ifstream infile(file);
  if (!infile.is_open()) {
    Rcpp::stop("failed to open input file");
  }
  auto dconv = create_distance_converter(threshold_config);
  Rcpp::IntegerMatrix im(dconv->m, seqnames.size());
  auto algo = create_cluster_algorithm(method_config, dconv.get())->create(im);
  auto worker = create_cluster_worker(parallel_config, algo.get(), infile);

  if (worker->n_threads() == 1) {
    (*worker)(0, 1);
  } else {
    RcppParallel::parallelFor(0, worker->n_threads(), *worker, 1, worker->n_threads());
  }
  worker->finalize();
  Rcpp::RObject output = R_NilValue;
  if (output_type == "matrix") {
    RcppParallel::RMatrix<int> im2(im);
    algo->write_to_matrix(im2);
    output = im;
  } else if (output_type == "hclust") {
    output = algo->as_hclust(seqnames);
  }
  return output;
}

// [[Rcpp::export]]
Rcpp::RObject distmx_cluster_multi(
    const std::string file,
    const Rcpp::CharacterVector seqnames,
    const Rcpp::ListOf<Rcpp::CharacterVector> which,
    const Rcpp::List threshold_config,
    const Rcpp::List method_config,
    const Rcpp::List parallel_config,
    const std::string output_type = "matrix"
) {
  if (output_type != "matrix" && output_type != "hclust") {
    Rcpp::stop("Unknown 'output_type'");
  }
  Rcpp::RObject output = R_NilValue;
  std::ifstream infile(file);
  if (!infile.is_open()) {
    Rcpp::stop("failed to open input file");
  }
  // Rcpp::Rcout << "creating DistanceConverter..." << std::flush;
  auto dconv = create_distance_converter(threshold_config);
  // Rcpp::Rcout << "done" << std::endl
              // << "creating ClusterAlgorithmFactory..." << std::flush;
  auto factory = create_cluster_algorithm(method_config, dconv.get());
  // Rcpp::Rcout << "done" << std::endl
  //             << "creating MultipleClusterAlgorithm..." << std::flush;
  auto algo = create_multiple_cluster_algorithm(parallel_config, *factory, seqnames, which);
  // Rcpp::Rcout << "done" << std::endl
  //             << "creating ClusterWorker..." << std::flush;
  auto worker = create_cluster_worker(parallel_config, algo.get(), infile);
  // Rcpp::Rcout << "done" << std::endl
  //             << "clustering..." << std::flush;

  if (worker->n_threads() == 1) {
    (*worker)(0, 1);
  } else {
    RcppParallel::parallelFor(0, worker->n_threads(), *worker, 1, worker->n_threads());
  }
  // Rcpp::Rcout << "done" << std::endl
  //             << "finalizing worker..." << std::flush;
  worker->finalize();
  // Rcpp::Rcout << "done" << std::endl
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
  return output;
}
