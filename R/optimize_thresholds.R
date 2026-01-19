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
  checkmate::assert_matrix(k, mode = "integer", row.names = "unique", any.missing = FALSE)
  checkmate::assert_numeric(as.integer(row.names(k)), any.missing = FALSE)
  thresholds <- as.numeric(row.names(k))
  checkmate::assert_integer(c, any.missing = FALSE)
  checkmate::assert_count(threads, positive = TRUE)
  checkmate::assert_character(measures)
  checkmate::assert_subset(measures, c("MCC", "RI", "ARI", "FMI", "MI", "AMI", "FM"))

  # measures which are based on the confusion matrix
  confusion_measures <- list(
    MCC = matthews_correlation_coefficient,
    RI = rand_index,
    ARI = adjusted_rand_index,
    FMI = fowlkes_mallow_index
  )
  confusion_measures <- confusion_measures[intersect(measures, names(confusion_measures))]

  # information theoretic measures
  information_measures <- intersect(c("MI", "AMI"), measures)

  # other measures
  other_measures <- list(FM = fmeasure)
  other_measures <- other_measures[intersect(names(other_measures), measures)]

  # the measures in the order they will be calculated
  sorted_measures <- c(names(confusion_measures), information_measures, names(other_measures))

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
  checkmate::assert_matrix(k, mode = "integer", row.names = "unique", any.missing = FALSE)
  checkmate::assert_numeric(as.integer(row.names(k)), any.missing = FALSE)
  thresholds <- as.numeric(row.names(k))
  checkmate::assert_integer(c, any.missing = FALSE)
  checkmate::assert_count(threads, positive = TRUE)
  checkmate::assert_character(measures)
  checkmate::assert_subset(measures, c("MCC", "RI", "ARI", "FMI", "MI", "AMI", "FM"))

  # measures which are based on the confusion matrix
  confusion_measures <- list(
    MCC = matthews_correlation_coefficient,
    RI = rand_index,
    ARI = adjusted_rand_index,
    FMI = fowlkes_mallow_index
  )
  confusion_measures <- confusion_measures[intersect(measures, names(confusion_measures))]

  # information theoretic measures
  information_measures <- intersect(c("MI", "AMI"), measures)

  # other measures
  other_measures <- list(FM = fmeasure)
  other_measures <- other_measures[intersect(names(other_measures), measures)]

  # the measures in the order they will be calculated
  sorted_measures <- c(names(confusion_measures), information_measures, names(other_measures))

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
#' [`DNAStringSet`][Biostrings::XStringSet-class], or `data.frame`) the
#' reference sequences
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
    ranks = c("kingdom", "phylum", "class", "order", "family", "genus",
               "species"),
    dist_config = dist_wfa2(),
    threshold_config = threshold_uniform(0.0, 0.4, 0.001),
    clust_config = clust_tree(),
    parallel_config = parallel_concurrent(threads = 1),
    min_taxa = 5L,
    min_refseq = 2L * min_taxa,
    id_col = "seq_id",
    measures = c("MCC", "RI", "ARI", "FMI", "MI", "AMI", "FM"),
    verbose = FALSE
) {

  # Check input
  checkmate::assert_character(ranks)
  checkmate::assert_string(id_col)
  checkmate::assert_data_frame(taxonomy)
  checkmate::assert_names(names(taxonomy), must.include = c(id_col, ranks))
  refseq_names <- seq_names(refseq)
  checkmate::assert_set_equal(taxonomy[[id_col]], refseq_names)
  checkmate::assert_class(threshold_config, "optimotu_threshold_config")
  checkmate::assert_class(dist_config, "optimotu_dist_config")
  checkmate::assert_class(parallel_config, "optimotu_parallel_config")
  checkmate::assert_integerish(min_refseq, lower = 1L)
  checkmate::assert_integerish(min_taxa, lower = 1L)
  checkmate::assert_character(measures)
  checkmate::assert_subset(measures, c("MCC", "RI", "ARI", "FMI", "MI", "AMI", "FM"))
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
  testset_select <- summarize_by_rank(taxonomy, ranks, id_col)
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

  # calculate clustering quality measures
  out <- mapply(
    find_best_threshold,
    c = testset_select$true_partition,
    k = clust,
    MoreArgs = list(
      threads = parallel_config$threads,
      measures = measures
    ),
    SIMPLIFY = FALSE
  )
  out <- do.call(rbind, out)
  cbind(
    rank = rep(testset_select$rank, each = length(measures)),
    superrank = rep(testset_select$superrank, each = length(measures)),
    supertaxon = rep(testset_select$supertaxon, each = length(measures)),
    out
  )
}
