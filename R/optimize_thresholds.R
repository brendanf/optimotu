# SPDX-CopyrightText: (c) 2025 Brendan Furneaux
# SPDX-License-Identifier: MIT

#' Identify "placeholder" taxa
#'
#' @param taxon (`character` vector) the taxon names to check
#' @return a logical vector indicating whether each taxon is a placeholder
#' @export
is_placeholder <- function(taxon) {
  checkmate::assert_character(taxon)
  is.na(taxon) |
    grepl(
      "^(dummy|unclassified|unknown|uncultured|environmental|unassigned|none)",
      taxon,
      ignore.case = TRUE
    ) |
    grepl("[_ ][Ss]p(\\b|_|[A-Z0-9])", taxon) |
    grepl("incertae[_ ]sedis", taxon, ignore.case = TRUE)
}

#' Replace "placeholder" ranks in a taxonomy with NA
#' @param taxonomy (`data.frame`) the taxonomy to clean
#' @param ranks (`character` vector) the ranks to clean
#' @return a cleaned taxonomy
#' @export
clean_taxonomy <- function(
  taxonomy,
  ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
) {
  checkmate::assert_data_frame(taxonomy)
  checkmate::assert_character(ranks)
  checkmate::assert_names(names(taxonomy), must.include = c("seq_id", ranks))

  for (rank in ranks) {
    taxonomy[[rank]] <-
      ifelse(is_placeholder(taxonomy[[rank]]), NA, taxonomy[[rank]])
  }
  taxonomy
}

#' Return the median of a (presorted) vector
#'
#' Always returns a value which is contained in `x`, even if `x` has an even
#' number of elements, in which case the lower of the two middle values is
#' returned.
#'
#' @param x (`numeric` vector) sorted vector to calculate the median of
#' @return the median of `x`
#' @export
#' @keywords internal
strict_median <- function(x) {
  if (length(x) == 0) {
    return(NA)
  } else if (length(x) == 1) {
    return(x)
  } else {
    return(x[length(x) %/% 2L + length(x) %% 2L])
  }
}

is_memory_budget_error <- function(e) {
  grepl("clustering memory budget exceeded", conditionMessage(e), fixed = TRUE)
}

estimate_subset_memory_mb <- function(
  n_seq,
  n_thresholds,
  clust_method,
  parallel_method,
  threads,
  include_result = FALSE
) {
  # Match Cluster*Factory::estimate_bytes() (int = 4, uint64_t = 8).
  int_bytes <- 4
  u64_bytes <- 8
  algo_bytes <- switch(
    clust_method,
    tree = n_seq * u64_bytes * 32,
    slink = n_seq * u64_bytes * 24,
    matrix = n_seq * n_thresholds * int_bytes + n_seq * int_bytes * 8,
    index = n_seq * n_thresholds * int_bytes + n_seq * int_bytes * 12,
    n_seq * n_thresholds * int_bytes
  )
  method_scale <- switch(
    parallel_method,
    merge = {
      n_tiles <- ceiling(0.5 * (sqrt(1 + 8 * threads) - 1))
      1 + threads * (2 / n_tiles)
    },
    concurrent = 1.1,
    hierarchical = max(1, threads / 2),
    1
  )
  total_bytes <- algo_bytes * method_scale
  # Tree/SLINK do not keep the n x m output; matrix/index already include it.
  if (isTRUE(include_result) && clust_method %in% c("tree", "slink")) {
    total_bytes <- total_bytes + n_seq * n_thresholds * int_bytes
  }
  max(1, total_bytes / (1024 * 1024))
}

estimate_batch_memory_mb <- function(
  batch_idx,
  testset_select,
  threshold_config,
  clust_config,
  parallel_config,
  include_result = FALSE
) {
  n_thresholds <- if (threshold_config$type == "uniform") {
    floor((threshold_config$to - threshold_config$from) / threshold_config$by) +
      1L
  } else {
    length(threshold_config$thresholds)
  }
  clust_method <- clust_config$method
  parallel_method <- parallel_config$method
  threads <- parallel_config$threads

  per_subset <- vapply(
    testset_select$n_seq[batch_idx],
    estimate_subset_memory_mb,
    numeric(1),
    n_thresholds = n_thresholds,
    clust_method = clust_method,
    parallel_method = parallel_method,
    threads = threads,
    include_result = include_result
  )
  sum(per_subset)
}

make_initial_batches <- function(
  ordered_idx,
  testset_select,
  threshold_config,
  clust_config,
  parallel_config,
  budget_mb
) {
  if (length(ordered_idx) == 0L) {
    return(list())
  }
  if (is.null(budget_mb) || !is.finite(budget_mb) || budget_mb <= 0) {
    return(list(ordered_idx))
  }

  batches <- list()
  current <- integer(0)
  for (idx in ordered_idx) {
    candidate <- c(current, idx)
    est <- estimate_batch_memory_mb(
      candidate,
      testset_select,
      threshold_config,
      clust_config,
      parallel_config
    )
    if (length(current) > 0L && est > budget_mb) {
      batches[[length(batches) + 1L]] <- current
      current <- idx
    } else {
      current <- candidate
    }
  }
  batches[[length(batches) + 1L]] <- current
  batches
}

split_batch_indices <- function(
  batch_idx,
  testset_select,
  strategy = "overlap"
) {
  if (length(batch_idx) <= 1L) {
    return(list(left = batch_idx, right = integer(0)))
  }
  if (identical(strategy, "overlap")) {
    groups <- split(
      batch_idx,
      testset_select$supertaxon[batch_idx],
      drop = TRUE
    )
    if (length(groups) > 1L) {
      group_sizes <- vapply(
        groups,
        function(g) sum(testset_select$n_seq[g]),
        numeric(1)
      )
      group_order <- names(sort(group_sizes, decreasing = TRUE))
      left <- integer(0)
      right <- integer(0)
      left_size <- 0
      right_size <- 0
      for (gname in group_order) {
        grp <- groups[[gname]]
        grp_size <- sum(testset_select$n_seq[grp])
        if (left_size <= right_size) {
          left <- c(left, grp)
          left_size <- left_size + grp_size
        } else {
          right <- c(right, grp)
          right_size <- right_size + grp_size
        }
      }
      if (length(left) > 0L && length(right) > 0L) {
        return(list(left = left, right = right))
      }
    }
  }
  o <- order(testset_select$n_seq[batch_idx], decreasing = TRUE, batch_idx)
  sorted <- batch_idx[o]
  mid <- length(sorted) %/% 2L
  list(left = sorted[seq_len(mid)], right = sorted[-seq_len(mid)])
}

set_parallel_threads <- function(parallel_config, threads) {
  out <- parallel_config
  out$threads <- as.integer(threads)
  out
}

set_parallel_subproblems <- function(parallel_config, subproblems) {
  out <- parallel_config
  out$subproblems <- as.integer(subproblems)
  out
}

#' Calculate clustering quality measures
#'
#' This function is primarily intended for use in plotting the clustering
#' quality measures at different thresholds. For choosing the best threshold,
#' use `optimize_thresholds()` or `find_best_threshold()`.
#'
#' The measures are abbreviated as follows:
#'
#'  - MCC: `matthews_correlation_coefficient()`
#'  - RI: `rand_index()`
#'  - ARI: `adjusted_rand_index()`
#'  - FMI: `fowlkes_mallow_index()`
#'  - MI: `adjusted_mutual_information()` (mutual information)
#'  - AMI: `adjusted_mutual_information()` (adjusted mutual information)
#'  - FM: `fmeasure()`
#'
#' A single call to this function with several measures may be faster than
#' multiple calls with a single measure, because the confusion matrix and/or
#' mutual information are calculated only once.
#'
#' @param k (`integer` matrix) the clustering partitions at different thresholds
#' @param c (`integer` vector) the true partition
#' @param threads (`integer(1)`) the number of threads to use
#' @param measures (`character` vector) the clustering quality measures to
#' calculate. Supported measures are "MCC", "RI", "ARI", "FMI", "MI", "AMI", and
#' "FM".
#' @return (`data.frame`) a data frame with the following columns:
#'  - `threshold` (`numeric`) the clustering threshold
#'  - `measure` (`character`) the clustering quality measure
#'  - `value` (`numeric`) the value of the measure at the threshold
#' The results are organized by the `measure`, and within each measure by the
#' thresholds in the same order as `k`.  However the `measure` column is not
#' guaranteed to be in the same order as the input.
#' @export
calculate_cluster_measures <- function(
  k,
  c,
  threads = 1L,
  measures = c("MCC", "RI", "ARI", "FMI", "MI", "AMI", "FM")
) {
  checkmate::assert_matrix(
    k,
    mode = "integer",
    row.names = "unique",
    any.missing = FALSE
  )
  checkmate::assert_numeric(as.integer(row.names(k)), any.missing = FALSE)
  thresholds <- as.numeric(row.names(k))
  checkmate::assert_integer(c, any.missing = FALSE)
  checkmate::assert_count(threads, positive = TRUE)
  checkmate::assert_character(measures)
  checkmate::assert_subset(
    measures,
    c("MCC", "RI", "ARI", "FMI", "MI", "AMI", "FM")
  )

  # measures which are based on the confusion matrix
  confusion_measures <- list(
    MCC = matthews_correlation_coefficient,
    RI = rand_index,
    ARI = adjusted_rand_index,
    FMI = fowlkes_mallow_index
  )
  confusion_measures <- confusion_measures[intersect(
    measures,
    names(confusion_measures)
  )]

  # information theoretic measures
  information_measures <- intersect(c("MI", "AMI"), measures)

  # other measures
  other_measures <- list(FM = fmeasure)
  other_measures <- other_measures[intersect(names(other_measures), measures)]

  # the measures in the order they will be calculated
  sorted_measures <- c(
    names(confusion_measures),
    information_measures,
    names(other_measures)
  )

  # allocate output
  n <- nrow(k) * length(measures)
  out <- data.frame(
    threshold = rep(thresholds, length(measures)),
    measure = rep(sorted_measures, each = nrow(k)),
    value = numeric(n)
  )

  # calculate values
  j <- 0L
  if (length(confusion_measures) > 0) {
    conf_mat <- confusion_matrix(k, c, threads)
    for (measure in names(confusion_measures)) {
      values <- confusion_measures[[measure]](conf_mat)
      out$value[j + seq_len(nrow(k))] <- values
      j <- j + nrow(k)
    }
  }
  if (length(information_measures) > 0) {
    mi <- adjusted_mutual_information(k, c, threads)
    for (measure in information_measures) {
      values <- mi[[measure]]
      out$value[j + seq_len(nrow(k))] <- values
      j <- j + nrow(k)
    }
  }
  for (measure in names(other_measures)) {
    out$value[j + seq_len(nrow(k))] <- other_measures[[measure]](k, c, threads)
    j <- j + nrow(k)
  }
  out
}

#' Find the optimal clustering threshold using one or more measures
#'
#' @param k (`integer` matrix) the clustering partitions at different thresholds
#' @param c (`integer` vector) the true partition
#' @param threads (`integer(1)`) the number of threads to use
#' @param measures (`character` vector) the clustering quality measures to
#' calculate. Supported measures are "MCC", "RI", "ARI", "FMI", "MI", "AMI", and
#' "FM".
#' @return (`data.frame`) a data frame with the following columns:
#' - `measure` (`character`) the clustering quality measure
#' - `threshold` (`numeric`) the threshold which maximizes the measure
#' - `value` (`numeric`) the value of the measure at the threshold
#' @export
find_best_threshold <- function(
  k,
  c,
  threads = 1L,
  measures = c("MCC", "RI", "ARI", "FMI", "MI", "AMI", "FM")
) {
  checkmate::assert_matrix(
    k,
    mode = "integer",
    row.names = "unique",
    any.missing = FALSE
  )
  checkmate::assert_numeric(as.integer(row.names(k)), any.missing = FALSE)
  thresholds <- as.numeric(row.names(k))
  checkmate::assert_integer(c, any.missing = FALSE)
  checkmate::assert_count(threads, positive = TRUE)
  checkmate::assert_character(measures)
  checkmate::assert_subset(
    measures,
    c("MCC", "RI", "ARI", "FMI", "MI", "AMI", "FM")
  )

  # measures which are based on the confusion matrix
  confusion_measures <- list(
    MCC = matthews_correlation_coefficient,
    RI = rand_index,
    ARI = adjusted_rand_index,
    FMI = fowlkes_mallow_index
  )
  confusion_measures <- confusion_measures[intersect(
    measures,
    names(confusion_measures)
  )]

  # information theoretic measures
  information_measures <- intersect(c("MI", "AMI"), measures)

  # other measures
  other_measures <- list(FM = fmeasure)
  other_measures <- other_measures[intersect(names(other_measures), measures)]

  # the measures in the order they will be calculated
  sorted_measures <- c(
    names(confusion_measures),
    information_measures,
    names(other_measures)
  )

  # allocate output
  out <- data.frame(
    measure = sorted_measures,
    threshold = NA_real_,
    value = NA_real_
  )

  # calculate values
  j <- 1L
  if (length(confusion_measures) > 0) {
    conf_mat <- confusion_matrix(k, c, threads)
    for (measure in names(confusion_measures)) {
      values <- confusion_measures[[measure]](conf_mat)
      best_value <- max(values, na.rm = TRUE)
      best_threshold <- thresholds[strict_median(which(values == best_value))]
      out$threshold[j] <- best_threshold
      out$value[j] <- best_value
      j <- j + 1L
    }
  }
  if (length(information_measures) > 0) {
    mi <- adjusted_mutual_information(k, c, threads)
    for (measure in information_measures) {
      values <- mi[[measure]]
      best_value <- max(values, na.rm = TRUE)
      best_threshold <- thresholds[strict_median(which(values == best_value))]
      out$threshold[j] <- best_threshold
      out$value[j] <- best_value
      j <- j + 1L
    }
  }
  for (measure in names(other_measures)) {
    values <- other_measures[[measure]](k, c, threads)
    best_value <- max(values, na.rm = TRUE)
    best_threshold <- thresholds[strict_median(which(values == best_value))]
    out$threshold[j] <- best_threshold
    out$value[j] <- best_value
    j <- j + 1L
  }
  out
}

#' Optimize clustering thresholds for taxonomically identified reference sequences
#'
#' @param taxonomy (`data.frame`) taxonomic identifications of the reference
#' sequences.  Must include a column for each rank in `ranks` and a column for
#' sequence identifiers, defined by `id_col`. Any additional columns are ignored.
#' @param refseq (named `character`, file name,
#' [`DNAStringSet`][Biostrings::XStringSet-class], `data.frame`,
#' `fastqindexr_index`, or character vector of `.fqi` paths) the
#' reference sequences
#' @param seq_idx optional 1-based subset indices for `refseq`; see [seq_as_char()].
#' @param seq_file optional path overrides for index-backed `refseq` only; see
#'   [seq_as_char()]. Taxonomy rows must match the selected sequence IDs.
#' @param ranks (`character` vector) the taxonomic ranks in `taxonomy`
#' @param dist_config (`op timotu_dist_config`) specification of the pairwise
#' distance algorithm to use, as created by `dist_config()` or its helpers
#' @param threshold_config (`optimotu_threshold_config`) specification of the
#' thresholds to test, as created by `threshold_config()` or its helpers
#' @param clust_config (`optimotu_clust_config`) specification of the clustering
#' algorithm to use, as created by `clust_config()` or its helpers
#' @param parallel_config (`optimotu_parallel_config`) specification of the
#' parallelization scheme to use, as created by `parallel_config()` or its
#' helpers
#' @param min_taxa (`integer(1)`) the minimum number of subtaxa at a given rank
#' which must belong to a taxon to optimize thresholds for that rank within that
#' taxon. Must be at least 2.
#' @param min_refseq (`integer(1)`) the minimum number of reference sequences
#' which must belong to a taxon to optimize thresholds for that taxon. Must be
#' at least `min_taxa` (but should probably be more).
#' @param id_col (`character(1)`) the name of the column in `taxonomy` which
#' contains sequence identifiers
#' @param measures (`character`) one or more measures to calculate optimum
#' thresholds for
#' @param clustering_memory_budget_mb (`numeric(1)` or `NULL`) Best-effort
#' memory budget for clustering-owned native structures, in MB. This is not a
#' hard process-memory cap; set below machine limits for safer behavior.
#' @param retry_on_memory_exhaustion (`logical(1)`) Retry clustering batches if
#' the clustering memory budget is exceeded.
#' @param retry_split_strategy (`character(1)`) Strategy used to split failed
#' batches. `"overlap"` prefers keeping related subsets together and
#' `"balanced"` splits by size.
#' @param retry_reduce_threads (`logical(1)`) Whether to reduce thread count
#' before splitting on budget failure.
#' @param retry_min_threads (`integer(1)`) Lower bound for thread-reduction
#' retries.
#' @param retry_overpartition (`logical(1)`) Whether to retry with more pair
#' subproblems than threads to reduce per-subproblem mapped state.
#' @param retry_subproblem_factor (`integer(1)`) Multiplicative growth factor
#' for subproblem count during overpartition retries.
#' @param verbose (`logical(1)` or `integer(1)`) whether to print progress
#' messages; values greater than 1 (or TRUE) print more
#' @return (`data.frame`) a data frame with the following columns:
#'   - `rank` (`character`) the rank being optimized
#'   - `superrank` (`character`) the containing rank for optimization
#'   - `supertaxon` (`character`) the containing taxon
#'   - `measure` (`character`) the clustering quality measure
#'   - `threshold` (`numeric`) the threshold which maximizes the measure
#'   - `value` (`numeric`) the value of the measure at the threshold
#' @export
optimize_thresholds <- function(
  taxonomy,
  refseq,
  ranks = c(
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  ),
  dist_config = dist_wfa2(),
  threshold_config = threshold_uniform(0.0, 0.4, 0.001),
  clust_config = clust_tree(),
  parallel_config = parallel_concurrent(threads = 1),
  min_taxa = 5L,
  min_refseq = 2L * min_taxa,
  id_col = "seq_id",
  measures = c("MCC", "RI", "ARI", "FMI", "MI", "AMI", "FM"),
  verbose = FALSE,
  seq_idx = NULL,
  seq_file = NULL,
  clustering_memory_budget_mb = NULL,
  retry_on_memory_exhaustion = TRUE,
  retry_split_strategy = c("overlap", "balanced"),
  retry_reduce_threads = TRUE,
  retry_min_threads = 1L,
  retry_overpartition = TRUE,
  retry_subproblem_factor = 2L
) {
  # Check input
  checkmate::assert_character(ranks)
  checkmate::assert_string(id_col)
  checkmate::assert_data_frame(taxonomy)
  checkmate::assert_names(names(taxonomy), must.include = c(id_col, ranks))
  refseq_names <- seq_input_seq_ids(refseq, seq_idx, seq_file)
  checkmate::assert_set_equal(taxonomy[[id_col]], refseq_names)
  checkmate::assert_class(threshold_config, "optimotu_threshold_config")
  checkmate::assert_class(dist_config, "optimotu_dist_config")
  checkmate::assert_class(parallel_config, "optimotu_parallel_config")
  checkmate::assert_integerish(min_refseq, lower = 1L)
  checkmate::assert_integerish(min_taxa, lower = 1L)
  checkmate::assert_character(measures)
  checkmate::assert_subset(
    measures,
    c("MCC", "RI", "ARI", "FMI", "MI", "AMI", "FM")
  )
  checkmate::assert_number(
    clustering_memory_budget_mb,
    lower = 1,
    null.ok = TRUE,
    finite = TRUE
  )
  checkmate::assert_flag(retry_on_memory_exhaustion)
  retry_split_strategy <- match.arg(retry_split_strategy)
  checkmate::assert_flag(retry_reduce_threads)
  checkmate::assert_count(retry_min_threads, positive = TRUE)
  checkmate::assert_flag(retry_overpartition)
  checkmate::assert_integerish(retry_subproblem_factor, lower = 2L, len = 1L)
  checkmate::assert(
    checkmate::check_flag(verbose),
    checkmate::check_integerish(verbose, lower = 0L)
  )

  # SLINK will fail if the order of sequences is different in the testset and
  # the reference sequences
  if (!isTRUE(all.equal(refseq_names, taxonomy[[id_col]]))) {
    taxonomy <- taxonomy[match(refseq_names, taxonomy[[id_col]]), ]
  }

  # Calculate which subsets to optimize
  if (isTRUE(verbose) || verbose >= 1L) {
    cat("Calculating subsets to optimize...")
  }
  testset_select <- summarize_by_rank(taxonomy, ranks, id_col)
  testset_select <- testset_select[
    (testset_select$n_seq >= min_refseq & testset_select$n_taxa >= min_taxa) |
      testset_select$supertaxon %in% taxonomy[[ranks[1]]],
  ]
  if (isTRUE(verbose) || verbose >= 1L) {
    cat(
      " done\n  Found",
      nrow(testset_select),
      "total thresholds to optimize.\n"
    )
    if (verbose >= 2L) {
      rank <- NULL
      superrank <- NULL
      dplyr::count(testset_select, rank, superrank) |>
        dplyr::mutate(
          dplyr::across(
            c(rank, superrank),
            \(x) rank2factor(x, ranks)
          )
        ) |>
        dplyr::arrange(rank, superrank) |>
        glue::glue_data(
          "  - {rank} within {n} {superrank}-rank taxa\n",
          .trim = FALSE
        ) |>
        cat(sep = "")
    }
  }

  row_order <- order(
    testset_select$rank,
    testset_select$superrank,
    testset_select$supertaxon,
    testset_select$n_seq,
    decreasing = TRUE
  )
  ordered_idx <- seq_len(nrow(testset_select))[row_order]

  initial_batches <- make_initial_batches(
    ordered_idx,
    testset_select,
    threshold_config,
    clust_config,
    parallel_config,
    clustering_memory_budget_mb
  )
  if (
    (isTRUE(verbose) || verbose >= 1L) && !is.null(clustering_memory_budget_mb)
  ) {
    cat(
      "Planning memory-aware batches:",
      length(initial_batches),
      "batch(es) for",
      length(ordered_idx),
      "subset problems.\n"
    )
  }

  run_batch <- function(batch_idx, local_parallel_config, depth = 0L) {
    if (length(batch_idx) == 0L) {
      return(NULL)
    }
    if (isTRUE(verbose) || verbose >= 2L) {
      cat(
        sprintf(
          "Running batch depth=%d n_subsets=%d threads=%d%s\n",
          depth,
          length(batch_idx),
          local_parallel_config$threads,
          if (!is.null(local_parallel_config$subproblems)) {
            paste0(" subproblems=", local_parallel_config$subproblems)
          } else {
            ""
          }
        )
      )
    }
    result <- tryCatch(
      {
        clust_args <- list(
          seq = refseq,
          dist_config = dist_config,
          threshold_config = threshold_config,
          clust_config = clust_config,
          parallel_config = local_parallel_config,
          output_type = "matrix",
          which = testset_select$seq_id[batch_idx],
          verbose = verbose,
          seq_idx = seq_idx,
          seq_file = seq_file
        )
        if (!is.null(clustering_memory_budget_mb)) {
          clust_args$clustering_memory_budget_mb <- clustering_memory_budget_mb
        }
        clust <- do.call(seq_cluster, clust_args)
        out <- mapply(
          find_best_threshold,
          c = testset_select$true_partition[batch_idx],
          k = clust,
          MoreArgs = list(
            threads = local_parallel_config$threads,
            measures = measures
          ),
          SIMPLIFY = FALSE
        )
        out <- do.call(rbind, out)
        cbind(
          row_index = rep(batch_idx, each = length(measures)),
          rank = rep(testset_select$rank[batch_idx], each = length(measures)),
          superrank = rep(
            testset_select$superrank[batch_idx],
            each = length(measures)
          ),
          supertaxon = rep(
            testset_select$supertaxon[batch_idx],
            each = length(measures)
          ),
          out
        )
      },
      error = identity
    )
    if (!inherits(result, "error")) {
      return(result)
    }
    if (!retry_on_memory_exhaustion || !is_memory_budget_error(result)) {
      stop(result)
    }
    if (
      retry_reduce_threads && local_parallel_config$threads > retry_min_threads
    ) {
      next_threads <- max(
        retry_min_threads,
        local_parallel_config$threads %/% 2L
      )
      if (next_threads < local_parallel_config$threads) {
        if (isTRUE(verbose) || verbose >= 1L) {
          cat(
            "Retrying batch with reduced threads:",
            local_parallel_config$threads,
            "->",
            next_threads,
            "\n"
          )
        }
        next_parallel <- set_parallel_threads(
          local_parallel_config,
          next_threads
        )
        return(run_batch(batch_idx, next_parallel, depth))
      }
    }
    if (
      retry_overpartition &&
        local_parallel_config$method %in% c("merge", "concurrent")
    ) {
      current_subproblems <- if (is.null(local_parallel_config$subproblems)) {
        local_parallel_config$threads
      } else {
        local_parallel_config$subproblems
      }
      next_subproblems <- max(
        current_subproblems * retry_subproblem_factor,
        length(batch_idx)
      )
      if (next_subproblems > current_subproblems) {
        if (isTRUE(verbose) || verbose >= 1L) {
          cat(
            "Retrying batch with more subproblems:",
            current_subproblems,
            "->",
            next_subproblems,
            "\n"
          )
        }
        next_parallel <- set_parallel_subproblems(
          local_parallel_config,
          next_subproblems
        )
        retry_result <- tryCatch(
          run_batch(batch_idx, next_parallel, depth),
          error = identity
        )
        if (!inherits(retry_result, "error")) {
          return(retry_result)
        }
        result <- retry_result
      }
    }
    if (length(batch_idx) == 1L) {
      stop(
        "Single subset still exceeds memory budget. Increase ",
        "`clustering_memory_budget_mb` or use lower-memory settings."
      )
    }
    split <- split_batch_indices(
      batch_idx,
      testset_select,
      strategy = retry_split_strategy
    )
    if (length(split$left) == 0L || length(split$right) == 0L) {
      split <- split_batch_indices(
        batch_idx,
        testset_select,
        strategy = "balanced"
      )
    }
    left_size <- sum(testset_select$n_seq[split$left])
    right_size <- sum(testset_select$n_seq[split$right])
    first <- if (left_size >= right_size) split$left else split$right
    second <- if (left_size >= right_size) split$right else split$left
    if (isTRUE(verbose) || verbose >= 1L) {
      cat(
        "Splitting batch at depth",
        depth,
        "into sizes",
        length(first),
        "and",
        length(second),
        "\n"
      )
    }
    out_first <- run_batch(first, local_parallel_config, depth + 1L)
    out_second <- run_batch(second, local_parallel_config, depth + 1L)
    rbind(out_first, out_second)
  }

  outputs <- lapply(
    initial_batches,
    run_batch,
    local_parallel_config = parallel_config
  )
  out <- do.call(rbind, outputs)
  out <- out[order(out$row_index), , drop = FALSE]
  out$row_index <- NULL
  out
}
