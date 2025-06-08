#ifdef OPTIMOTU_R

#include "config.h"

//' Search for best match(es) of query sequences in reference sequences
//' @param query (named `character`) query sequences
//' @param ref (named `character`) reference sequences
//' @param dist_config (`optimotu_dist_config`) distance configuration
//' @param parallel_config (`optimotu_parallel_config`) parallelization
//' configuration. Search methods do not treat the different methods differently;
//' only the number of threads is used.
//' @param threshold (`numeric`) distance threshold for the search; larger
//' distances are not considered
//' @param verbose (`integer`) verbosity level
//' @return `data.frame` with columns `seq_id` (query sequence ID),
//' `ref_id` (reference sequence ID), and `dist` (distance between the query
//' and reference sequences)
//' @export
//' @keywords internal
// [[Rcpp::export]]
Rcpp::RObject seq_search_internal(
  Rcpp::CharacterVector query,
  Rcpp::CharacterVector ref,
  Rcpp::List dist_config,
  Rcpp::List parallel_config,
  double threshold,
  int verbose = 0,
  bool return_cigar = false,
  int span = 0
) {
  if (!query.hasAttribute("names") || !ref.hasAttribute("names")) {
    Rcpp::stop("query and ref must be named");
  }
  Rcpp::CharacterVector query_names = query.names();
  Rcpp::CharacterVector ref_names = ref.names();

  auto query_str = Rcpp::as<std::vector<std::string>>(query);
  auto ref_str = Rcpp::as<std::vector<std::string>>(ref);
  auto search_worker = create_search_worker(
    dist_config, parallel_config, query_str, ref_str, threshold, verbose,
    return_cigar, span
  );
  int threads = search_worker->n_threads();
  if (threads == 1) {
    (*search_worker)(0, 1);
  } else {
    RcppParallel::parallelFor(0, threads, *search_worker, 1, threads);
  }

  // count size of output
  int n = 0;
  for (auto & hit : search_worker->get_hits()) {
    n += hit->best_ref.size();
  }

  // allocate outputs
  Rcpp::CharacterVector query_out(n);
  Rcpp::CharacterVector ref_out(n);
  Rcpp::NumericVector dist_out(n);
  Rcpp::CharacterVector cigar_out(n);

  // fill outputs
  int i = 0, j = 0;
  for (auto & hits : search_worker->get_hits()) {
    int start_i = i;
    for (auto & hit : hits->best_ref) {
      query_out[i] = query_names[j];
      ref_out[i] = ref_names[hit];
      dist_out[i] = hits->best_dist;
      ++i;
    }
    if (return_cigar) {
      i = start_i;
      for (auto & cigar : hits->best_cigar()) {
        cigar_out[i] = cigar;
        ++i;
      }
    }
    ++j;
  }

  // create output
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("seq_id") = query_out,
    Rcpp::Named("ref_id") = ref_out,
    Rcpp::Named("dist") = dist_out
  );
  if (return_cigar) {
    out["cigar"] = cigar_out;
  }
  out.attr("class") = Rcpp::CharacterVector::create("tbl_df", "tbl", "data.frame");
  out.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -n);
  return out;
}


#endif
