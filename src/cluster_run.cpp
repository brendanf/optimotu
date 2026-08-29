// SPDX-FileCopyrightText: 2025 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifdef OPTIMOTU_R

#include "cluster_run.h"
#include "optimotu.h"
#include <RcppParallel.h>

void stop_memory_budget(const MemoryBudgetExceeded &e) {
  const double budget_mb =
    static_cast<double>(e.budget_bytes) / (1024.0 * 1024.0);
  const double current_mb =
    static_cast<double>(e.current_bytes) / (1024.0 * 1024.0);
  const double requested_mb =
    static_cast<double>(e.requested_bytes) / (1024.0 * 1024.0);
  OPTIMOTU_STOP(
    "clustering memory budget exceeded (budget=%.2f MB, current=%.2f MB, "
    "requested=%.2f MB, context=%s)",
    budget_mb, current_mb, requested_mb, e.context.c_str()
  );
}

MultiClusterJob run_seq_cluster_multi(
  const SequenceSet &seq,
  Rcpp::CharacterVector seqnames,
  Rcpp::ListOf<Rcpp::CharacterVector> which,
  Rcpp::List dist_config,
  Rcpp::List threshold_config,
  Rcpp::List clust_config,
  Rcpp::List parallel_config,
  int verbose,
  double clustering_memory_budget_mb
) {
  MultiClusterJob job;
  if (verbose) {
    OPTIMOTU_CERR << "creating DistanceConverter..." << std::flush;
  }
  job.dconv = create_distance_converter(threshold_config);
  if (verbose) {
    OPTIMOTU_CERR << "done\ncreating ClusterAlgorithmFactory..."
                  << std::flush;
  }
  job.factory = create_cluster_algorithm(clust_config, job.dconv.get());
  if (verbose) {
    OPTIMOTU_CERR << "done\ncreating MultipleClusterAlgorithm..."
                  << std::flush;
  }
  std::size_t clustering_memory_budget_bytes = 0;
  if (clustering_memory_budget_mb > 0.0) {
    clustering_memory_budget_bytes = static_cast<std::size_t>(
      clustering_memory_budget_mb * 1024.0 * 1024.0
    );
  }
  job.algo = create_multiple_cluster_algorithm(
    parallel_config,
    *job.factory,
    seqnames,
    which,
    verbose,
    clustering_memory_budget_bytes
  );
  if (verbose) {
    OPTIMOTU_CERR << "done\ncreating ClusterWorker..." << std::flush;
  }
  try {
    if (seq.size() >= 2) {
      auto worker = create_dist_cluster_worker(
        dist_config, parallel_config, seq, *job.algo, verbose
      );
      if (verbose) {
        OPTIMOTU_CERR << "done\nclustering..." << std::endl;
      }
      if (worker->n_threads() == 1) {
        (*worker)(0, 1);
      } else {
        RcppParallel::parallelFor(
          0, worker->n_threads(), *worker, 1, worker->n_threads()
        );
      }
    } else if (verbose) {
      OPTIMOTU_CERR << "done\nskipping ClusterWorker for " << seq.size()
                    << " input sequence(s)" << std::endl;
    }
  } catch (const MemoryBudgetExceeded &e) {
    stop_memory_budget(e);
  }
  if (verbose) {
    OPTIMOTU_CERR << "done\nfinalizing worker..." << std::flush;
  }
  job.algo->finalize();
  job.algo->prepare_output();
  return job;
}

MultiClusterJob run_seq_cluster_multi(
  const SequenceSet &seq,
  Rcpp::ListOf<Rcpp::IntegerVector> which,
  Rcpp::List dist_config,
  Rcpp::List threshold_config,
  Rcpp::List clust_config,
  Rcpp::List parallel_config,
  int verbose,
  double clustering_memory_budget_mb
) {
  MultiClusterJob job;
  if (verbose) {
    OPTIMOTU_CERR << "creating DistanceConverter..." << std::flush;
  }
  job.dconv = create_distance_converter(threshold_config);
  if (verbose) {
    OPTIMOTU_CERR << "done\ncreating ClusterAlgorithmFactory..."
                  << std::flush;
  }
  job.factory = create_cluster_algorithm(clust_config, job.dconv.get());
  if (verbose) {
    OPTIMOTU_CERR << "done\ncreating MultipleClusterAlgorithm..."
                  << std::flush;
  }
  std::size_t clustering_memory_budget_bytes = 0;
  if (clustering_memory_budget_mb > 0.0) {
    clustering_memory_budget_bytes = static_cast<std::size_t>(
      clustering_memory_budget_mb * 1024.0 * 1024.0
    );
  }
  job.algo = create_multiple_cluster_algorithm(
    parallel_config,
    *job.factory,
    static_cast<j_t>(seq.size()),
    which,
    verbose,
    clustering_memory_budget_bytes
  );
  if (verbose) {
    OPTIMOTU_CERR << "done\ncreating ClusterWorker..." << std::flush;
  }
  try {
    if (seq.size() >= 2) {
      auto worker = create_dist_cluster_worker(
        dist_config, parallel_config, seq, *job.algo, verbose
      );
      if (verbose) {
        OPTIMOTU_CERR << "done\nclustering..." << std::endl;
      }
      if (worker->n_threads() == 1) {
        (*worker)(0, 1);
      } else {
        RcppParallel::parallelFor(
          0, worker->n_threads(), *worker, 1, worker->n_threads()
        );
      }
    } else if (verbose) {
      OPTIMOTU_CERR << "done\nskipping ClusterWorker for " << seq.size()
                    << " input sequence(s)" << std::endl;
    }
  } catch (const MemoryBudgetExceeded &e) {
    stop_memory_budget(e);
  }
  if (verbose) {
    OPTIMOTU_CERR << "done\nfinalizing worker..." << std::flush;
  }
  job.algo->finalize();
  job.algo->prepare_output();
  return job;
}

MultiClusterJob run_distmx_cluster_multi(
  std::istream &infile,
  Rcpp::CharacterVector seqnames,
  Rcpp::ListOf<Rcpp::CharacterVector> which,
  Rcpp::List threshold_config,
  Rcpp::List clust_config,
  Rcpp::List parallel_config,
  bool by_name,
  bool verbose
) {
  MultiClusterJob job;
  if (verbose) {
    OPTIMOTU_CERR << "creating DistanceConverter..." << std::flush;
  }
  job.dconv = create_distance_converter(threshold_config);
  if (verbose) {
    OPTIMOTU_CERR << "done\ncreating ClusterAlgorithmFactory..."
                  << std::flush;
  }
  job.factory = create_cluster_algorithm(clust_config, job.dconv.get());
  if (verbose) {
    OPTIMOTU_CERR << "done\ncreating MultipleClusterAlgorithm..."
                  << std::flush;
  }
  job.algo = create_multiple_cluster_algorithm(
    parallel_config, *job.factory, seqnames, which
  );
  if (verbose) {
    OPTIMOTU_CERR << "done\ncreating ClusterWorker..." << std::flush;
  }
  auto worker = create_cluster_worker<std::istream>(
    parallel_config, job.algo.get(), infile, by_name, seqnames
  );
  if (verbose) {
    OPTIMOTU_CERR << "done\nclustering..." << std::flush;
  }
  if (worker->n_threads() == 1) {
    (*worker)(0, 1);
  } else {
    RcppParallel::parallelFor(
      0, worker->n_threads(), *worker, 1, worker->n_threads()
    );
  }
  if (verbose) {
    OPTIMOTU_CERR << "done\nfinalizing worker..." << std::flush;
  }
  worker->finalize();
  job.algo->finalize();
  job.algo->prepare_output();
  return job;
}

#endif // OPTIMOTU_R
