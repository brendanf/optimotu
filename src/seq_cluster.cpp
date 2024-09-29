#ifdef OPTIMOTU_R

#include "config.h"

// [[Rcpp::export]]
Rcpp::RObject seq_cluster_single(
    const Rcpp::CharacterVector &seq,
    const Rcpp::List dist_config,
    const Rcpp::List threshold_config,
    const Rcpp::List clust_config,
    const Rcpp::List parallel_config,
    const std::string output_type = "matrix",
    const bool verbose = false
) {
  if (output_type != "matrix" && output_type != "hclust") {
    OPTIMOTU_STOP("Unknown 'output_type'");
  }
  auto dconv = create_distance_converter(threshold_config);
  Rcpp::IntegerMatrix im(dconv->m, seq.size());
  std::vector<std::string> cppseq = Rcpp::as<std::vector<std::string>>(seq);
  auto algo = create_cluster_algorithm(clust_config, dconv.get())->create(im);
  auto worker = create_align_cluster_worker(dist_config, parallel_config, cppseq, *algo, verbose);
  int threads = worker->n_threads();
  if (threads == 1) {
    (*worker)(0, 1);
  } else {
    RcppParallel::parallelFor(0, threads, *worker, 1, threads);
  }
  if (verbose) {
    OPTIMOTU_COUT << "back in main function" << std::endl;
    OPTIMOTU_COUT << worker->aligned() << " aligned / "
                  << worker->prealigned() << " prealigned"
                  << std::endl;
  }
  algo->finalize();
  Rcpp::RObject output = R_NilValue;
  if (output_type == "matrix") {
    auto method = (Rcpp::as<Rcpp::CharacterVector>(clust_config["method"]))[0];
    if (method == "tree" || method == "slink") {
      internal_matrix_t m(im);
      algo->write_to_matrix(m);
    }
    output = im;
  } else if (output_type == "hclust") {
    output = algo->as_hclust(seq.names());
  }
  return output;
}

// [[Rcpp::export]]
Rcpp::RObject seq_cluster_multi(
    const Rcpp::CharacterVector &seq,
    const Rcpp::ListOf<Rcpp::CharacterVector> which,
    const Rcpp::List dist_config,
    const Rcpp::List threshold_config,
    const Rcpp::List clust_config,
    const Rcpp::List parallel_config,
    const std::string output_type = "matrix",
    const bool verbose = false
) {
  if (output_type != "matrix" && output_type != "hclust") {
    OPTIMOTU_STOP("Unknown 'output_type'");
  }
  Rcpp::RObject output = R_NilValue;

  const std::vector<std::string> cppseq = Rcpp::as<std::vector<std::string>>(seq);

  // OPTIMOTU_COUT << "creating DistanceConverter..." << std::flush;
  auto dconv = create_distance_converter(threshold_config);
  // OPTIMOTU_COUT << "done" << std::endl
  // << "creating ClusterAlgorithmFactory..." << std::flush;
  auto factory = create_cluster_algorithm(clust_config, dconv.get());
  // OPTIMOTU_COUT << "done" << std::endl
  //             << "creating MultipleClusterAlgorithm..." << std::flush;
  auto algo = create_multiple_cluster_algorithm(parallel_config, *factory, seq.names(), which);
  // OPTIMOTU_COUT << "done" << std::endl
  //             << "creating ClusterWorker..." << std::flush;
  auto worker = create_align_cluster_worker(dist_config, parallel_config, cppseq, *algo, verbose);
  // OPTIMOTU_COUT << "done" << std::endl
  //             << "clustering..." << std::flush;

  if (worker->n_threads() == 1) {
    (*worker)(0, 1);
  } else {
    RcppParallel::parallelFor(0, worker->n_threads(), *worker, 1, worker->n_threads());
  }
  // OPTIMOTU_COUT << "done" << std::endl
  //             << "finalizing worker..." << std::flush;
  algo->finalize();
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
  return output;
}
#endif // OPTIMOTU_R
