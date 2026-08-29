#ifdef OPTIMOTU_R

#include "config.h"
#include "cluster_run.h"
#include "cluster_measures.h"
#include <RcppParallel.h>
#include <algorithm>
#include <chrono>
#include <numeric>
#include <random>
#include <vector>

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
  size_t n_prealigned = 0, n_aligned = 0;
  if (cppseq.size() >= 2) {
    if (verbose) {
      OPTIMOTU_CERR << "done\ncreating ClusterWorker..." << std::flush;
    }
    auto worker = create_dist_cluster_worker(dist_config, parallel_config, cppseq, *algo, verbose);
    if (verbose) {
      OPTIMOTU_CERR << "done\nclustering..." << std::endl;
    }
    int threads = worker->n_threads();
    if (threads == 1) {
      (*worker)(0, 1);
    } else {
      RcppParallel::parallelFor(0, threads, *worker, 1, threads);
    }
    n_prealigned = worker->prealigned();
    n_aligned = worker->aligned();
  } else if (verbose) {
    OPTIMOTU_CERR << "done\nskipping ClusterWorker for " << cppseq.size()
                  << " input sequence(s)\n";
  }
  if (verbose) {
    OPTIMOTU_CERR << "done\n"
                  << n_aligned << " aligned / "
                  << n_prealigned << " prealigned"
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
    const int verbose = 0,
    const double clustering_memory_budget_mb = -1.0
) {
  if (output_type != "matrix" && output_type != "hclust") {
    OPTIMOTU_STOP("Unknown 'output_type'");
  }
  const std::vector<std::string> cppseq = Rcpp::as<std::vector<std::string>>(seq);
  auto job = run_seq_cluster_multi(
      cppseq,
      seq.names(),
      which,
      dist_config,
      threshold_config,
      clust_config,
      parallel_config,
      verbose,
      clustering_memory_budget_mb);
  if (verbose)
    OPTIMOTU_CERR << "done\ncreating output..." << std::flush;
  Rcpp::List output;
  if (output_type == "matrix") {
    Rcpp::NumericVector thresh(job.dconv->m);
    for (int i = 0; i < job.dconv->m; ++i)
    {
      thresh[i] = job.dconv->inverse(i);
    }

    auto outlist = Rcpp::List(which.size());
    auto internal_out = std::vector<RcppParallel::RMatrix<int>>();
    internal_out.reserve(which.size());
    for (int i = 0; i < which.size(); ++i) {
      auto outi = Rcpp::IntegerMatrix(job.dconv->m, which[i].size());
      outi.attr("dimnames") = Rcpp::List::create(
        thresh,
        which[i]
      );
      outlist[i] = outi;
      internal_out.emplace_back(outi);
    }
    job.algo->write_to_matrix(internal_out);
    output = outlist;
  } else if (output_type == "hclust") {
    output = job.algo->as_hclust();
  }
  output.names() = which.names();
  if (verbose)
    OPTIMOTU_CERR << "done" << std::endl;
  return output;
}

// [[Rcpp::export]]
Rcpp::List seq_cluster_multi_via_rows(
    const Rcpp::CharacterVector &seq,
    const Rcpp::ListOf<Rcpp::CharacterVector> which,
    const Rcpp::List dist_config,
    const Rcpp::List threshold_config,
    const Rcpp::List clust_config,
    const Rcpp::List parallel_config,
    const int verbose = 0,
    const double clustering_memory_budget_mb = -1.0)
{
  const std::vector<std::string> cppseq = Rcpp::as<std::vector<std::string>>(seq);
  auto job = run_seq_cluster_multi(
      cppseq,
      seq.names(),
      which,
      dist_config,
      threshold_config,
      clust_config,
      parallel_config,
      verbose,
      clustering_memory_budget_mb);
  Rcpp::NumericVector thresh(job.dconv->m);
  for (int i = 0; i < job.dconv->m; ++i)
  {
    thresh[i] = job.dconv->inverse(i);
  }
  Rcpp::List output(which.size());
  for (int i = 0; i < which.size(); ++i)
  {
    auto outi = subset_matrix_via_rows(*job.algo, static_cast<std::size_t>(i));
    outi.attr("dimnames") = Rcpp::List::create(thresh, which[i]);
    output[i] = outi;
  }
  output.names() = which.names();
  return output;
}

// [[Rcpp::export]]
Rcpp::List seq_cluster_multi_best_threshold(
    const Rcpp::CharacterVector &seq,
    const Rcpp::ListOf<Rcpp::CharacterVector> which,
    const Rcpp::List dist_config,
    const Rcpp::List threshold_config,
    const Rcpp::List clust_config,
    const Rcpp::List parallel_config,
    const Rcpp::ListOf<Rcpp::IntegerVector> true_partitions,
    const std::vector<std::string> &measures,
    const Rcpp::NumericVector &thresholds,
    const Rcpp::IntegerVector &threshold_order,
    const int verbose = 0,
    const double clustering_memory_budget_mb = -1.0)
{
  if (true_partitions.size() != which.size())
  {
    OPTIMOTU_STOP("true_partitions and which must have the same length");
  }
  const std::vector<std::string> cppseq = Rcpp::as<std::vector<std::string>>(seq);
  auto job = run_seq_cluster_multi(
      cppseq,
      seq.names(),
      which,
      dist_config,
      threshold_config,
      clust_config,
      parallel_config,
      verbose,
      clustering_memory_budget_mb);
  int threads = element_as_int(parallel_config, "threads", "parallel_config");
  Rcpp::List output(which.size());
  for (int i = 0; i < which.size(); ++i)
  {
    output[i] = find_best_threshold_from_algo(
        *job.algo,
        static_cast<std::size_t>(i),
        true_partitions[i],
        measures,
        thresholds,
        threshold_order,
        threads);
  }
  output.names() = which.names();
  return output;
}

namespace
{

  double elapsed_seconds(
      std::chrono::steady_clock::time_point t0,
      std::chrono::steady_clock::time_point t1)
  {
    return std::chrono::duration<double>(t1 - t0).count();
  }

  struct ThresholdRowFillWorker : public RcppParallel::Worker
  {
    SingleClusterAlgorithm &algo;
    explicit ThresholdRowFillWorker(SingleClusterAlgorithm &algo) : algo(algo) {}
    void operator()(std::size_t begin, std::size_t end)
    {
      std::vector<int> buf(static_cast<std::size_t>(algo.n));
      for (std::size_t t = begin; t < end; ++t)
      {
        algo.write_threshold_row(static_cast<d_t>(t), buf.data());
      }
    }
  };

} // namespace

// Time clustering vs write_threshold_row / write_to_matrix after clustering.
// Unexported; used to check that random-access row fill is cheap.
// [[Rcpp::export]]
Rcpp::DataFrame seq_cluster_profile_output(
    const Rcpp::CharacterVector &seq,
    const Rcpp::List dist_config,
    const Rcpp::List threshold_config,
    const Rcpp::List clust_config,
    const Rcpp::List parallel_config,
    const int row_repeats = 5,
    const int verbose = 0)
{
  if (row_repeats < 1)
  {
    OPTIMOTU_STOP("row_repeats must be at least 1");
  }
  if (!clust_config.inherits("optimotu_cluster_config"))
  {
    OPTIMOTU_STOP("clust_config is not a valid cluster configuration");
  }
  if (!threshold_config.inherits("optimotu_threshold_config"))
  {
    OPTIMOTU_STOP("threshold_config is not a valid threshold configuration");
  }
  if (!parallel_config.inherits("optimotu_parallel_config"))
  {
    OPTIMOTU_STOP("parallel_config is not a valid parallel configuration");
  }
  auto dconv = create_distance_converter(threshold_config);
  std::vector<std::string> cppseq = Rcpp::as<std::vector<std::string>>(seq);
  const j_t n = static_cast<j_t>(cppseq.size());
  const d_t m = dconv->m;
  auto factory = create_cluster_algorithm(clust_config, dconv.get());
  auto algo = factory->create(n);
  int threads = element_as_int(parallel_config, "threads", "parallel_config");
  if (threads < 1)
    threads = 1;

  auto t_cluster0 = std::chrono::steady_clock::now();
  if (cppseq.size() >= 2)
  {
    auto worker = create_dist_cluster_worker(
        dist_config, parallel_config, cppseq, *algo, verbose);
    int wthreads = worker->n_threads();
    if (wthreads == 1)
    {
      (*worker)(0, 1);
    }
    else
    {
      RcppParallel::parallelFor(0, wthreads, *worker, 1, wthreads);
    }
  }
  auto t_cluster1 = std::chrono::steady_clock::now();

  auto t_final0 = std::chrono::steady_clock::now();
  algo->finalize();
  auto t_final1 = std::chrono::steady_clock::now();

  auto t_prep0 = std::chrono::steady_clock::now();
  algo->prepare_output();
  auto t_prep1 = std::chrono::steady_clock::now();

  std::vector<int> dest(static_cast<std::size_t>(n));
  std::vector<d_t> order(static_cast<std::size_t>(m));
  std::iota(order.begin(), order.end(), 0);
  std::mt19937 rng(1);
  std::shuffle(order.begin(), order.end(), rng);

  auto t_seq0 = std::chrono::steady_clock::now();
  for (int r = 0; r < row_repeats; ++r)
  {
    for (d_t t = 0; t < m; ++t)
    {
      algo->write_threshold_row(t, dest.data());
    }
  }
  auto t_seq1 = std::chrono::steady_clock::now();

  auto t_shuf0 = std::chrono::steady_clock::now();
  for (int r = 0; r < row_repeats; ++r)
  {
    for (d_t t : order)
    {
      algo->write_threshold_row(t, dest.data());
    }
  }
  auto t_shuf1 = std::chrono::steady_clock::now();

  ThresholdRowFillWorker row_worker(*algo);
  auto t_par0 = std::chrono::steady_clock::now();
  for (int r = 0; r < row_repeats; ++r)
  {
    if (threads <= 1)
    {
      row_worker(0, static_cast<std::size_t>(m));
    }
    else
    {
      RcppParallel::parallelFor(
          0, static_cast<std::size_t>(m), row_worker, 1, threads);
    }
  }
  auto t_par1 = std::chrono::steady_clock::now();

  Rcpp::IntegerMatrix im(m, n);
  internal_matrix_t mat(im);
  auto t_mat0 = std::chrono::steady_clock::now();
  algo->write_to_matrix(mat);
  auto t_mat1 = std::chrono::steady_clock::now();

  const double repeats = static_cast<double>(row_repeats);
  Rcpp::CharacterVector phase = Rcpp::CharacterVector::create(
      "clustering",
      "finalize",
      "prepare_output",
      "write_threshold_row_sequential",
      "write_threshold_row_shuffled",
      "write_threshold_row_parallel",
      "write_to_matrix");
  Rcpp::NumericVector seconds = Rcpp::NumericVector::create(
      elapsed_seconds(t_cluster0, t_cluster1),
      elapsed_seconds(t_final0, t_final1),
      elapsed_seconds(t_prep0, t_prep1),
      elapsed_seconds(t_seq0, t_seq1) / repeats,
      elapsed_seconds(t_shuf0, t_shuf1) / repeats,
      elapsed_seconds(t_par0, t_par1) / repeats,
      elapsed_seconds(t_mat0, t_mat1));
  Rcpp::IntegerVector n_repeats = Rcpp::IntegerVector::create(
      1, 1, 1, row_repeats, row_repeats, row_repeats, 1);
  Rcpp::DataFrame out = Rcpp::DataFrame::create(
      Rcpp::Named("phase") = phase,
      Rcpp::Named("seconds") = seconds,
      Rcpp::Named("repeats") = n_repeats);
  out.attr("n_seq") = static_cast<int>(n);
  out.attr("n_thresholds") = static_cast<int>(m);
  out.attr("threads") = threads;
  out.attr("checksum") = dest.empty() ? 0 : dest[0];
  return out;
}
#endif // OPTIMOTU_R
