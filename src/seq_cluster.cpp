#ifdef OPTIMOTU_R

#include "config.h"

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix seq_cluster(
    const Rcpp::CharacterVector &seq,
    const Rcpp::List dist_config,
    const Rcpp::List threshold_config,
    const Rcpp::List clust_config,
    const Rcpp::List parallel_config,
    const bool verbose = false
) {
  auto dconv = create_distance_converter(threshold_config);
  Rcpp::IntegerMatrix im(dconv->m, seq.size());
  std::vector<std::string> cppseq = Rcpp::as<std::vector<std::string>>(seq);
  auto cm = create_cluster_algorithm(clust_config, dconv.get())->create(im);
  auto worker = create_align_cluster_worker(dist_config, parallel_config, cppseq, *cm, verbose);
  int threads = worker->nthreads();
  if (threads > 1) {
    RcppParallel::parallelFor(0, threads, *worker, 1, threads);
  } else {
    (*worker)(0, 1);
  }
  if (verbose) {
    std::cout << "back in main function" << std::endl;
    std::cout << worker->aligned() << " aligned / "
              << worker->prealigned() << " prealigned"
              << std::endl;
  }
  cm->finalize();
  auto method = (Rcpp::as<Rcpp::CharacterVector>(clust_config["method"]))[0];
  if (method == "tree" || method == "slink") {
    internal_matrix_t m(im);
    cm->write_to_matrix(m);
  }
  return im;
}

#endif
