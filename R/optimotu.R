# SPDX-CopyrightText: (c) 2025, Brendan Furneaux
# SPDX-License-Identifier: MIT

#' Perform closed-reference pseudo-single-linkage clustering by iterative
#' searching
#'
#' When a sequence is equally distant from multiple references, a random
#' reference is chosen.
#'
#' @param query (file name, named `character`, `data.frame` or
#' [`DNAStringSet`][Biostrings::XStringSet-class]) sequences to search for
#' @param ref (file name, named `character`, `data.frame` or
#' [`DNAStringSet`][Biostrings::XStringSet-class]) sequences to search within
#' @param threshold (`numeric`) distance threshold for clustering
#' @param ... additional arguments to pass to `seq_search()`
#' @return (`data.frame`) table with three columns: `seq_id` for the ID of the
#' query sequence, `ref_id` for the ID of the reference sequence.
#' @export
closed_ref_cluster <- function(query, ref, threshold, ...) {
  last_out <- data.frame(seq_id = character(0), ref_id = character(0))
  out = list(last_out)
  while (sequence_size(ref) > 0 && sequence_size(query) > 0) {
    result <- seq_search(query = query, ref = ref, threshold = threshold, ...)
    if (nrow(result) > 0) {
      result <- result[result$dist <= threshold, c("seq_id", "ref_id"), drop = FALSE]
      if (nrow(result) > 0) {
        result <- split(result, result$seq_id)
        result <- lapply(result, function(x) {
          if (nrow(x) > 1) {
            x <- x[sample.int(nrow(x), 1), ]
          }
          x
        })
        result <- do.call(rbind, result)
      }
    }
    if (nrow(last_out) == 0) {
      last_out <- result <- result[, c("seq_id", "ref_id")]
    } else {
      result <- result[, c("seq_id", "ref_id")]
      names(result)[2] <- "temp"
      result <- merge(
        result,
        last_out,
        by.x = "temp",
        by.y = "seq_id"
      )
      last_out <- result <- result[, c("seq_id", "ref_id")]
    }
    out <- c(out, list(result))
    ref <- select_sequence(query, result$seq_id)
    query <- select_sequence(query, result$seq_id, negate = TRUE)
  }
  do.call(rbind, out)
}

#' Taxonomically guided OTU clustering using optimized thresholds
#'
#' @param seqs ([`DNAStringSet`][Biostrings::XStringSet-class], file name,
#' `character` vector, or `data.frame`) input sequences
#' @param tax_prob (`data.frame`) taxonomic probabilities for each sequence,
#' with columns `seq_id` giving unique sequence identifiers, `rank` giving the
#' taxonomic rank, `taxon` giving the identified taxon name, and, optionally,
#' `prob` giving the probability of the taxon assignment. All sequences in `seqs`
#' must be identified at the first rank listed in `ranks` (usually "kingdom");
#' if this is not the case, it is recommended to add an additional dummy rank,
#' such as "root_rank", to `tax_prob` with the appropriate taxon assignments.
#' Rows where `taxon` or `prob` is `NA`, or where `prob < prob_thresh`, are
#' considered unknown.
#' @param threshold_optima (`data.frame`) optimized thresholds, as calculated
#' by `optimize_thresholds()`
#' @param ranks (`character`) The taxonomic ranks to cluster at, in order from
#' most inclusive to least inclusive. Defaults to `c("kingdom", "phylum",
#' "class", "order", "family", "genus", "species")`. Any extra ranks given in
#' `tax_prob` or `threshold_optima` are silently ignored.
#' @param measure (`character`) The measure to optimize for. Required if
#' multiple measures are present in `threshold_optima`, in which case the given
#' measure must be included in `threshold_optima$measure`.
#' @param prob_thresh (`numeric`) The probability threshold to use as a cutoff
#' for accepting taxonomic assignments. Defaults to `NULL`, in which case the
#' all taxonomic assignments in `tax_prob` are used. If given, `tax_prob` must
#' include a `prob` column.
#' @param dist_config (`optimotu_dist_config`) configuration for calculating
#' distances, as returned by `dist_config()` or its helpers.
#' @param clust_config (`optimotu_clust_config`) configuration for clustering,
#' as returned by `clust_config()` or its helpers.
#' @param parallel_config (`optimotu_parallel_config`) configuration for parallel
#' processing, as returned by `parallel_config()` or its helpers.
#' @param verbose (`logical` or `integer`) print progress messages.
#' @return (`data.frame`) with columns `seq_id` giving the sequence identifiers
#' and the taxonomic assignments at each rank in `ranks`
#' @export
optimotu <- function(
    seqs,
    tax_prob,
    threshold_optima,
    ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
    measure = NULL,
    prob_thresh = NULL,
    dist_config = dist_wfa2(),
    clust_config = clust_slink(),
    parallel_config = parallel_concurrent(1),
    verbose = FALSE
) {
  # check inputs
  checkmate::assert(
    checkmate::test_character(seqs, any.missing = FALSE),
    checkmate::test_data_frame(seqs),
    checkmate::test_class(seqs, "XStringSet")
  )
  # convert seqs to a data.frame if it isn't already
  if (is.character(seqs)) {
    if (length(seqs) == 1 && file.exists(seqs)) {
      if (grepl(fasta_regex, seqs)) {
        seqs <- Biostrings::readDNAStringSet(seqs)
      } else if (grepl(fastq_regex, seqs)) {
        seqs <- Biostrings::readDNAStringSet(seqs, format = "fastq")
      } else {
        stop("Input file must be in FASTA or FASTQ format")
      }
      seqs <- data.frame(seq_id = names(seqs), seq = as.character(seqs))
    } else {
      checkmate::assert_named(seqs, type = "unique")
      seqs <- data.frame(seq_id = names(seqs), seq = as.character(seqs))
    }
  } else if (methods::is(seq, "XStringSet")) {
    checkmate::assert_named(seqs, type = "unique")
    seqs <- data.frame(seq_id = names(seqs), seq = as.character(seqs))
  }

  checkmate::assert_data_frame(seqs)
  checkmate::assert_names(names(seqs), must.include = c("seq_id", "seq"))
  checkmate::assert_character(seqs$seq_id, unique = TRUE, any.missing = FALSE)
  checkmate::assert_character(seqs$seq, any.missing = FALSE)

  checkmate::assert_data_frame(tax_prob)
  checkmate::assert_names(names(tax_prob), must.include = c("seq_id", "rank", "taxon"))
  checkmate::assert_character(tax_prob$seq_id, any.missing = FALSE)
  checkmate::assert_set_equal(tax_prob$seq_id, seqs$seq_id, )
  checkmate::assert_character(tax_prob$rank, any.missing = FALSE)
  checkmate::assert_character(tax_prob$taxon)

  checkmate::assert_data_frame(threshold_optima)
  checkmate::assert_names(
    names(threshold_optima),
    must.include = c("rank", "superrank", "supertaxon", "threshold")
  )

  checkmate::assert_character(ranks)
  checkmate::assert_string(measure, null.ok = TRUE)
  if (!is.null(measure)) {
    checkmate::assert_subset(measure, threshold_optima$measure)
  }

  checkmate::assert_number(prob_thresh, null.ok = TRUE)
  if (!is.null(prob_thresh)) {
    checkmate::assert_names(names(tax_prob), must.include = "prob")
    checkmate::assert_numeric(tax_prob$prob)
  }

  checkmate::assert_class(dist_config, "optimotu_dist_config")
  if (identical(dist_config$method, "file")) {
    if (isFALSE(dist_config$by_name)) {
      stop("File-based distance matrices indexed by integer are not supported ",
           "in optimotu(). If your file has sequence names, use ",
           "dist_file({your_file}, by_name = TRUE) instead.")
    }
  }

  checkmate::assert_class(parallel_config, "optimotu_parallel_config")

  checkmate::assert(
    checkmate::test_flag(verbose),
    checkmate::test_count(verbose)
  )
  verbose <- as.integer(verbose)

  n_ranks <- length(ranks)

  # initialize the loop with kingdom-level "results"
  # these are the only variables which are actually needed between iterations
  taxon_table <- tax_prob[tax_prob$rank == ranks[1], c("seq_id", "taxon")]
  names(taxon_table)[2] <- ranks[1]
  checkmate::assert_character(taxon_table[[ranks[1]]], any.missing = FALSE)

  pseudotaxon_table <- data.frame(seq_id = character(0))
  pseudotaxon_table[[ranks[1]]] <- character(0)

  for (r in seq_len(n_ranks - 1)) {
    .rank <- ranks[r + 1]
    .parent_rank <- ranks[r]

    # these will have been updated in the previous iteration
    .parent_taxa <- taxon_table
    .parent_pseudotaxa <- pseudotaxon_table

    # define known taxonomy table for the rank in question
    if (is.null(prob_thresh)) {
      known_taxa <- tax_prob$rank == .rank & !is.na(tax_prob$taxon)
    } else {
      known_taxa <- tax_prob$rank == .rank & !is.na(tax_prob$taxon) &
        tax_prob$prob >= prob_thresh
    }
    known_taxon_table <- tax_prob[known_taxa, c("seq_id", "taxon")]
    names(known_taxon_table)[2] <- .rank
    known_taxon_table <- merge(
      .parent_taxa,
      known_taxon_table,
      by = "seq_id",
      all.x = TRUE
    )

    # Find groups that need to be closed-reference clustered
    preclosed_taxon_table <- split(known_taxon_table, known_taxon_table[[.parent_rank]])
    to_closedref_cluster <- vapply(
      preclosed_taxon_table,
      function(x) {any(is.na(x[[.rank]])) & !all(is.na(x[[.rank]]))},
      logical(1)
    )
    preclosed_taxon_table <- preclosed_taxon_table[to_closedref_cluster]

    # define the thresholds for that rank based on the parent rank
    thresholds <- calc_taxon_thresholds(
      rank = .parent_rank,
      taxon_table = known_taxon_table,
      optima = threshold_optima,
      ranks = ranks,
      measure = measure
    )

    # closed-reference clustering
    clusters_closed_ref <-lapply(
        preclosed_taxon_table,
        function(parent_group) {
          if (verbose) {
            cat("Closed reference clustering within", .parent_rank,
                parent_group[[.parent_rank]], "\n")
          }
          unknowns <- is.na(parent_group[[.rank]])

          taxon <- parent_group[[.parent_rank]][1]

          clusters <- closed_ref_cluster(
            query = select_sequence(seqs, parent_group$seq_id[unknowns]),
            ref = select_sequence(seqs, parent_group$seq_id[!unknowns]),
            threshold = threshold_as_dist(thresholds[taxon]),
            dist_config = dist_config,
            parallel_config = parallel_config,
            verbose = max(verbose - 1L, 0L)
          )

          # take taxonomy from the reference sequence
          new_taxa <- merge(
            clusters,
            parent_group[, c("seq_id", .rank)],
            by.x = "ref_id",
            by.y = "seq_id"
          )

          new_taxa[, c("seq_id", .rank)]
        }
      )
    clusters_closed_ref <- do.call(rbind, clusters_closed_ref)

    # add the clustering results to the taxon table
    known_taxon_table[match(clusters_closed_ref$seq_id, known_taxon_table$seq_id), .rank] <-
      clusters_closed_ref[[.rank]]

    # make a taxonomy table for use in denovo clustering
    predenovo_taxon_table <- known_taxon_table[is.na(known_taxon_table[[.rank]]),]

    # calculate the thresholds for denovo clustering ####
    denovo_thresholds <- calc_subtaxon_thresholds(
      rank = .parent_rank,
      taxon_table = predenovo_taxon_table,
      optima = threshold_optima,
      ranks = ranks,
      measure = measure
    )

    # split the pre-denovo taxon table by the parent rank
    predenovo_taxon_table <- split(predenovo_taxon_table, predenovo_taxon_table[[.parent_rank]])

    # denovo clustering
    clusters_denovo <- lapply(
      predenovo_taxon_table,
      function(parent_group) {
        .parent_taxon <- parent_group[[.parent_rank]][1]
        if (verbose) {
          cat("De novo clustering within", .parent_rank, .parent_taxon, "\n")
        }
        if (nrow(parent_group) > 1) {
          clusters <- seq_cluster(
            seq = select_sequence(seqs, parent_group$seq_id),
            threshold_config = threshold_set(
              if (.parent_taxon %in% names(denovo_thresholds)) {
                denovo_thresholds[[.parent_taxon]]
              } else {
                denovo_thresholds[["_NA_"]]
              }
            ),
            dist_config = dist_config,
            clust_config = clust_config,
            parallel_config = parallel_config,
            output_type = "matrix",
            verbose = max(0L, verbose - 1L)
          )
          clusters <- t(clusters)
          clusters <- as.data.frame(clusters, stringsAsFactors = FALSE)
          for (new_rank in names(clusters)) {
            parent_group[[new_rank]] <- clusters[[new_rank]]
          }
        } else {
          for (new_rank in c(.rank, subranks(.rank, ranks))) {
            parent_group[[new_rank]] <- 0L
          }
        }
        parent_group
      }
    )
    clusters_denovo <- do.call(rbind, clusters_denovo)

    if (is.null (clusters_denovo) || nrow(clusters_denovo) == 0) {
      clusters_denovo <- data.frame(seq_id = character(0))
      for (new_rank in superranks(.rank, ranks)) {
        clusters_denovo[[new_rank]] <- character(0)
      }
      for (new_rank in c(.rank, subranks(.rank, ranks))) {
        clusters_denovo[[new_rank]] <- integer(0)
      }
    }

    #### Step: make a taxonomy table to all assigned ranks ####
    taxon_table <- known_taxon_table[!is.na(known_taxon_table[[.rank]]),,drop = FALSE]

    #### Step: add the pseudotaxon names to the de novo clustered ASVs at current rank ####
    pseudotaxon_table <- rbind(clusters_denovo, .parent_pseudotaxa)
    pseudotaxon_table <- sort_by(pseudotaxon_table, ~ seq_id)
    pseudotaxon_table[[.rank]] <- paste(
      pseudotaxon_table[[.parent_rank]],
      pseudotaxon_table[[.rank]]
    )
    pseudotaxon_table[[.rank]] <- ordered(
      pseudotaxon_table[[.rank]],
      levels = unique(pseudotaxon_table[[.rank]]),
      labels = make_seq_names(
        n = length(unique(pseudotaxon_table[[.rank]])),
        prefix = paste0("pseudo", .rank, "_")
      )
    ) |> as.character()
  }
  # export the outputs to a final table
  out <- sort_by(
    rbind(taxon_table, pseudotaxon_table),
    ~ seq_id
  )
  rownames(out) <- NULL
  class(out) <- c("tbl_df", "tbl", "data.frame")
  out
}
