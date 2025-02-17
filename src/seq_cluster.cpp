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
    const int verbose = 0
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
  if (verbose) {
    OPTIMOTU_CERR << "creating DistanceConverter..." << std::flush;
  }
  auto dconv = create_distance_converter(threshold_config);
  if (verbose) {
    OPTIMOTU_CERR << "done\ncreating ClusterAlgorithm..." << std::flush;
  }
  Rcpp::IntegerMatrix im(dconv->m, seq.size());
  std::vector<std::string> cppseq = Rcpp::as<std::vector<std::string>>(seq);
  auto algo = create_cluster_algorithm(clust_config, dconv.get())->create(im);
  if (verbose) {
    OPTIMOTU_CERR << "done\ncreating ClusterWorker..." << std::flush;
  }
  auto worker = create_align_cluster_worker(dist_config, parallel_config, cppseq, *algo, verbose);
  if (verbose) {
    OPTIMOTU_CERR << "done\nclustering..." << std::endl;
  }
  int threads = worker->n_threads();
  if (threads == 1) {
    (*worker)(0, 1);
  } else {
    RcppParallel::parallelFor(0, threads, *worker, 1, threads);
  }
  if (verbose) {
    OPTIMOTU_CERR << "done\n"
                  << worker->aligned() << " aligned / "
                  << worker->prealigned() << " prealigned"
                  << "\ncreating output..." << std::flush;
  }
  algo->finalize();
  Rcpp::RObject output = R_NilValue;
  if (output_type == "matrix") {
    auto method = element_as_string(clust_config, "method", "clust_config");
    if (method == "tree" || method == "slink") {
      internal_matrix_t m(im);
      algo->write_to_matrix(m);
    }
    output = im;
    Rcpp::NumericVector thresh(dconv->m);
    for (int i = 0; i < dconv->m; ++i) {
      thresh[i] = dconv->inverse(i);
    }

    output.attr("dimnames") = Rcpp::List::create(
      thresh,
      seq.names()
    );
  } else if (output_type == "hclust") {
    output = algo->as_hclust(seq.names());
  }
  if (verbose) {
    OPTIMOTU_CERR << "done" << std::endl;
  }
  return output;
}

// [[Rcpp::export]]
Rcpp::List seq_cluster_multi(
    const Rcpp::CharacterVector &seq,
    const Rcpp::ListOf<Rcpp::CharacterVector> which,
    const Rcpp::List dist_config,
    const Rcpp::List threshold_config,
    const Rcpp::List clust_config,
    const Rcpp::List parallel_config,
    const std::string output_type = "matrix",
    const int verbose = 0
) {
  if (output_type != "matrix" && output_type != "hclust") {
    OPTIMOTU_STOP("Unknown 'output_type'");
  }
  Rcpp::List output;

  const std::vector<std::string> cppseq = Rcpp::as<std::vector<std::string>>(seq);

  if (verbose)
    OPTIMOTU_CERR << "creating DistanceConverter..." << std::flush;
  auto dconv = create_distance_converter(threshold_config);
  if (verbose) {
    OPTIMOTU_CERR << "done\ncreating ClusterAlgorithmFactory..." << std::flush;
  }
  auto factory = create_cluster_algorithm(clust_config, dconv.get());
  if (verbose) {
    OPTIMOTU_CERR << "done\ncreating MultipleClusterAlgorithm..." << std::flush;
  }
  auto algo = create_multiple_cluster_algorithm(parallel_config, *factory, seq.names(), which);
  if (verbose)
    OPTIMOTU_CERR << "done\ncreating ClusterWorker..." << std::flush;
  auto worker = create_align_cluster_worker(dist_config, parallel_config, cppseq, *algo, verbose);
  if (verbose)
    OPTIMOTU_CERR << "done\nclustering..." << std::endl;

  if (worker->n_threads() == 1) {
    (*worker)(0, 1);
  } else {
    RcppParallel::parallelFor(0, worker->n_threads(), *worker, 1, worker->n_threads());
  }
  if (verbose)
    OPTIMOTU_CERR << "done\nfinalizing worker..." << std::flush;
  algo->finalize();
  if (verbose)
    OPTIMOTU_CERR << "done\ncreating output..." << std::flush;
  if (output_type == "matrix") {
    Rcpp::NumericVector thresh(dconv->m);
    for (int i = 0; i < dconv->m; ++i) {
      thresh[i] = dconv->inverse(i);
    }

    auto outlist = Rcpp::List(which.size());
    auto internal_out = std::vector<RcppParallel::RMatrix<int>>();
    internal_out.reserve(which.size());
    for (int i = 0; i < which.size(); ++i) {
      auto outi = Rcpp::IntegerMatrix(dconv->m, which[i].size());
      outi.attr("dimnames") = Rcpp::List::create(
        thresh,
        which[i]
      );
      outlist[i] = outi;
      internal_out.emplace_back(outi);
    }
    algo->write_to_matrix(internal_out);
    output = outlist;
  } else if (output_type == "hclust") {
    output = algo->as_hclust();
  }
  output.names() = which.names();
  if (verbose)
    OPTIMOTU_CERR << "done" << std::endl;
  return output;
}
#endif // OPTIMOTU_R
