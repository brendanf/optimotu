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
  grepl("^(dummy|unclassified|unknown|uncultured|environmental|unassigned|none)", taxon, ignore.case = TRUE) |
    grepl("[_ ][Ss]p(\\b|_|[A-Z0-9])", taxon) |
    grepl("incertae[_ ]sedis", taxon, ignore.case = TRUE)
}

#' Replace "placeholder" ranks in a taxonomy with NA
#' @param taxonomy (`data.frame`) the taxonomy to clean
#' @param ranks (`character` vector) the ranks to clean
#' @return a cleaned taxonomy
#' @export
clean_taxonomy <- function(taxonomy, ranks) {
  checkmate::assert_data_frame(taxonomy)
  checkmate::assert_character(ranks)
  checkmate::assert_names(names(taxonomy), must.include = c("seq_id", ranks))

  for (rank in ranks) {
    taxonomy[[rank]] <-
      ifelse(is_placeholder(taxonomy[[rank]]), NA, taxonomy[[rank]])
  }
  taxonomy
}

#' Optimize clustering thresholds for taxonomically identified reference sequences
#'
#' @param taxonomy (`data.frame`) taxonomic identifications of the reference
#' sequences.  Must include a column for each rank in `ranks` and a column for
#' sequence identifiers, defined by `id_col`. Any additional columns are ignored.
#' @param refseq (named `character`) the reference sequences
#' @param ranks (`character` vector) the taxonomic ranks in `taxonomy`
#' @param dist_config (`optimotu_dist_config`) specification of the pairwise
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
#' @param verbose (`logical(1)` or `integer(1)`) whether to print progress
#' messages; values greater than 1 (or TRUE) print more
#' @return (`data.frame`) a data frame with the following columns:
#'   - `rank` (`character`) the rank being optimized
#'   - `superrank` (`character`) the containing rank for optimization
#'   - `supertaxon` (`character`) the containing taxon
#'   - `metric` (`character`) the clustering quality metric
#'   - `threshold` (`numeric`) the threshold which maximizes the metric
#'   - `value` (`numeric`) the value of the metric at the threshold
#' @export
optimize_thresholds <- function(
    taxonomy,
    refseq,
    ranks = c("kingdom", "phylum", "class", "order", "family", "genus",
               "species"),
    dist_config = dist_wfa2(),
    threshold_config = threshold_uniform(0.0, 0.4, 0.001),
    clust_config = clust_tree(),
    parallel_config = parallel_concurrent(threads = 1),
    min_taxa = 5L,
    min_refseq = 2L * min_taxa,
    id_col = "seq_id",
    verbose = FALSE
) {

  # Check input
  checkmate::assert_character(ranks)
  checkmate::assert_string(id_col)
  checkmate::assert_data_frame(taxonomy)
  checkmate::assert_names(names(taxonomy), must.include = c(id_col, ranks))
  checkmate::assert_character(refseq)
  checkmate::assert_set_equal(taxonomy[[id_col]], names(refseq))
  checkmate::assert_class(threshold_config, "optimotu_threshold_config")
  checkmate::assert_class(dist_config, "optimotu_dist_config")
  checkmate::assert_class(parallel_config, "optimotu_parallel_config")
  checkmate::assert_integerish(min_refseq, lower = 1L)
  checkmate::assert_integerish(min_taxa, lower = 1L)
  checkmate::assert(
    checkmate::check_flag(verbose),
    checkmate::check_integerish(verbose, lower = 0L)
  )

  # SLINK will fail if the order of sequences is different in the testset and
  # the reference sequences
  if (!isTRUE(all.equal(names(refseq), taxonomy[[id_col]]))) {
    taxonomy <- taxonomy[match(names(refseq), taxonomy[[id_col]]), ]
  }

  # Calculate which subsets to optimize
  testset_select <- summarize_by_rank(taxonomy, ranks)
  testset_select <- testset_select[
    (testset_select$n_seq >= min_refseq & testset_select$n_taxa >= min_taxa) |
      testset_select$supertaxon %in% taxonomy[[ranks[1]]],
  ]

  # Do test clustering
  clust <- seq_cluster(
    refseq,
    dist_config = dist_config,
    threshold_config = threshold_config,
    clust_config = clust_config,
    parallel_config = parallel_config,
    output_type = "matrix",
    which = testset_select$seq_id,
    verbose = verbose
  )

  n_metrics <- 7L

  # calculate clustering quality measures
  out <- data.frame(
    rank = character(nrow(testset_select) * n_metrics),
    superrank = NA_character_,
    supertaxon = NA_character_,
    metric = NA_character_,
    threshold = NA_real_,
    value = NA_real_
  )
  j <- 1L
  for (i in seq_len(nrow(testset_select))) {
    conf_mat <- confusion_matrix(k = clust[[i]], c = testset_select$true_taxa[[i]], threads = parallel_config$threads)
    confusion_metrics <- list(
      MCC = matthews_correlation_coefficient,
      RI = rand_index,
      ARI = adjusted_rand_index,
      FMI = fowlkes_mallow_index
    )
    for (metric in names(confusion_metrics)) {
      value <- confusion_metrics[[metric]](conf_mat)
      best_value <- max(value, na.rm = TRUE)
      best_threshold <- rownames(clust[[i]])[utils::tail(n = 1, which(value == best_value))]
      out$rank[j] <- testset_select$rank[[i]]
      out$superrank[j] <- testset_select$superrank[[i]]
      out$supertaxon[j] <- testset_select$supertaxon[[i]]
      out$metric[j] <- metric
      out$threshold[j] <- as.numeric(best_threshold)
      out$value[j] <- best_value
      j <- j + 1L
    }
    mi <- adjusted_mutual_information(
      k = clust[[i]],
      c = testset_select$true_taxa[[i]],
      threads = parallel_config$threads
    )
    information_metrics <- c("MI", "AMI")
    for (metric in information_metrics) {
      value <- mi[[metric]]
      best_value <- max(value, na.rm = TRUE)
      best_threshold <- rownames(clust[[i]])[utils::tail(n = 1, which(value == best_value))]
      out$rank[j] <- testset_select$rank[[i]]
      out$superrank[j] <- testset_select$superrank[[i]]
      out$supertaxon[j] <- testset_select$supertaxon[[i]]
      out$metric[j] <- metric
      out$threshold[j] <- as.numeric(best_threshold)
      out$value[j] <- best_value
      j <- j + 1L
    }
    fm <- fmeasure(
      k = clust[[i]],
      c = testset_select$true_taxa[[i]],
      ncpu = parallel_config$threads
    )
    best_value <- max(fm, na.rm = TRUE)
    best_threshold <- rownames(clust[[i]])[utils::tail(n = 1, which(fm == best_value))]
    out$rank[j] <- testset_select$rank[[i]]
    out$superrank[j] <- testset_select$superrank[[i]]
    out$supertaxon[j] <- testset_select$supertaxon[[i]]
    out$metric[j] <- "FM"
    out$threshold[j] <- as.numeric(best_threshold)
    out$value[j] <- best_value
    j <- j + 1L
  }
  out
}
