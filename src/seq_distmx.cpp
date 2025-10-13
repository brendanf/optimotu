#ifdef OPTIMOTU_R

#include "config.h"

//' Search for best match(es) of query sequences in reference sequences
//'
//' @param seq (named `character`) sequences
//' @param dist_config (`optimotu_dist_config`) distance configuration
//' @param parallel_config (`optimotu_parallel_config`) parallelization
//' configuration. Search methods do not treat the different methods differently;
//' only the number of threads is used.
//' @param threshold (`numeric`) distance threshold for the search; larger
//' distances are not considered
//' @param verbose (`integer`) verbosity level
//' @param details (`integer`) detail level; 0 = score only, 1 = score and gap
//' stats, 2 = score and cigar
//' @param span (`integer`) alignment span
//' @return `data.frame` with columns `seq_id1` (first sequence ID),
//' `seq_id2` (second sequence ID), `score1` and `score2` (optimal alignment
//' score for the two sequences in the "prealignment" (if any) and "alignment"
//' stages), `dist1` and `dist2` (corresponding sequence distances). If
//' `details` is 1, also includes columns `align_length`, `n_insert`,
//' `n_delete`, `max_insert`, and `max_delete`. If `details` is 2, also
//' includes column `cigar`.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::RObject seq_distmx_internal(
  Rcpp::CharacterVector seq,
  Rcpp::List dist_config,
  Rcpp::List parallel_config,
  double threshold,
  int verbose = 0,
  int details = 0,
  int span = 0,
  bool constrain = true
) {

  // allocate output vectors
  std::vector<size_t> seq_id1;
  std::vector<size_t> seq_id2;
  std::vector<int> score1;
  std::vector<int> score2;
  std::vector<double> dist1;
  std::vector<double> dist2;
  // outputs for gap stats
  std::vector<int> aln_len;
  std::vector<int> n_ins;
  std::vector<int> n_del;
  std::vector<int> max_ins;
  std::vector<int> max_del;
  // outputs for cigar
  std::vector<std::string> cigar;

  auto seq_str = Rcpp::as<std::vector<std::string>>(seq);
  std::unique_ptr<SparseDistanceMatrix> sdm;
  if (details == 0) {
    sdm = std::make_unique<SparseDistanceMatrix>(seq_id1, seq_id2, score1, score2, dist1, dist2);
  } else if (details == 1) {
    sdm = std::make_unique<SparseDistanceMatrixGapstats>(seq_id1, seq_id2, score1, score2, dist1, dist2, aln_len, n_ins, n_del, max_ins, max_del);
  } else if (details == 2) {
    sdm = std::make_unique<SparseDistanceMatrixCigar>(seq_id1, seq_id2, score1, score2, dist1, dist2, cigar);
  } else {
    OPTIMOTU_STOP("details must be 0, 1, or 2");
  }
  auto align_worker = create_dist_worker(
    dist_config,
    parallel_config,
    seq_str,
    threshold,
    *sdm,
    verbose,
    span,
    constrain
  );
  int threads = align_worker->n_threads();
  if (threads == 1) {
    (*align_worker)(0, 1);
  } else {
    RcppParallel::parallelFor(0, threads, *align_worker, 1, threads);
  }

  // count size of output
  int n = sdm->seq1.size();

  // allocate outputs
  Rcpp::IntegerVector seq1_out(n);
  Rcpp::IntegerVector seq2_out(n);
  Rcpp::NumericVector score1_out(n);
  Rcpp::NumericVector score2_out(n);
  Rcpp::NumericVector dist1_out(n);
  Rcpp::NumericVector dist2_out(n);

  // fill outputs
  for (int i = 0; i < n; i++) {
    seq1_out[i] = sdm->seq1[i] + 1;
    seq2_out[i] = sdm->seq2[i] + 1;
    score1_out[i] = sdm->score1[i];
    score2_out[i] = sdm->score2[i];
    dist1_out[i] = sdm->dist1[i];
    dist2_out[i] = sdm->dist2[i];
  }

  // create output
  Rcpp::List out = Rcpp::List::create(
    Rcpp::Named("seq_idx1") = seq1_out,
    Rcpp::Named("seq_idx2") = seq2_out,
    Rcpp::Named("score1") = score1_out,
    Rcpp::Named("score2") = score2_out,
    Rcpp::Named("dist1") = dist1_out,
    Rcpp::Named("dist2") = dist2_out
  );

  if (details == 1) {
    auto sdm_gapstats = dynamic_cast<SparseDistanceMatrixGapstats*>(sdm.get());
    Rcpp::IntegerVector align_length_out(n);
    Rcpp::IntegerVector n_insert_out(n);
    Rcpp::IntegerVector n_delete_out(n);
    Rcpp::IntegerVector max_insert_out(n);
    Rcpp::IntegerVector max_delete_out(n);
    for (int i = 0; i < n; i++) {
      align_length_out[i] = sdm_gapstats->align_length[i];
      n_insert_out[i] = sdm_gapstats->n_insert[i];
      n_delete_out[i] = sdm_gapstats->n_delete[i];
      max_insert_out[i] = sdm_gapstats->max_insert[i];
      max_delete_out[i] = sdm_gapstats->max_delete[i];
    }
    out["align_length"] = align_length_out;
    out["n_insert"] = n_insert_out;
    out["n_delete"] = n_delete_out;
    out["max_insert"] = max_insert_out;
    out["max_delete"] = max_delete_out;
  }

  if (details == 2) {
    auto sdm_cigar = dynamic_cast<SparseDistanceMatrixCigar*>(sdm.get());
    Rcpp::CharacterVector cigar_out(n);
    for (int i = 0; i < n; i++) {
      cigar_out[i] = sdm_cigar->cigar[i];
    }
    out["cigar"] = cigar_out;
  }

  if (verbose >= 1) {
    Rcpp::Rcerr << "Prealigned: " << align_worker->prealigned()
                << "\nAligned: " << align_worker->aligned() << std::endl;
  }

  // set attributes
  out.attr("class") = Rcpp::CharacterVector::create("tbl_df", "tbl", "data.frame");
  out.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -n);
  return out;
}

#endif
