testthat::test_that("taxonomy_sequence_order is lexical by ranks then id", {
  taxonomy <- data.frame(
    seq_id = c("s3", "s1", "s2", "s4"),
    kingdom = c("k1", "k1", "k1", "k1"),
    phylum = c("p2", "p1", "p1", "p2"),
    class = c("c3", "c1", "c2", "c4"),
    stringsAsFactors = FALSE
  )
  ord <- optimotu:::taxonomy_sequence_order(
    taxonomy,
    ranks = c("kingdom", "phylum", "class"),
    id_col = "seq_id"
  )
  sorted <- taxonomy[ord, , drop = FALSE]
  testthat::expect_equal(sorted$seq_id, c("s1", "s2", "s3", "s4"))
  # Members of phylum p1 are a contiguous block.
  p1 <- which(sorted$phylum == "p1")
  testthat::expect_equal(p1, seq.int(min(p1), max(p1)))
})

testthat::test_that("optimize_thresholds is invariant to caller sequence order", {
  seq <- c(
    "ACGTACGTACGT",
    "ACGTACGTTCGT",
    "ACGTTCGTACGT",
    "ACGTTCGTTCGT",
    "ACGTACGGACGT",
    "ACGTACGGTCGT"
  )
  names(seq) <- paste0("s", seq_along(seq))
  taxonomy <- data.frame(
    seq_id = names(seq),
    kingdom = rep("k1", length(seq)),
    phylum = rep(c("p1", "p1", "p1", "p2", "p2", "p2"), 1),
    class = rep(c("c1", "c1", "c2", "c3", "c3", "c4"), 1),
    stringsAsFactors = FALSE
  )

  # Already in taxonomic order.
  sorted_tax <- taxonomy[
    optimotu:::taxonomy_sequence_order(
      taxonomy,
      ranks = c("kingdom", "phylum", "class"),
      id_col = "seq_id"
    ),
    ,
    drop = FALSE
  ]
  sorted_seq <- seq[sorted_tax$seq_id]

  set.seed(1L)
  shuf <- sample(seq_along(seq))
  shuffled_seq <- seq[shuf]
  shuffled_tax <- taxonomy[shuf, , drop = FALSE]

  common_args <- list(
    ranks = c("kingdom", "phylum", "class"),
    dist_config = dist_hamming(),
    threshold_config = threshold_set(c(0.05, 0.10, 0.20)),
    clust_config = clust_slink(),
    parallel_config = parallel_merge(threads = 2L),
    min_taxa = 2L,
    min_refseq = 3L,
    measures = c("MCC", "FM"),
    verbose = FALSE
  )

  from_sorted <- do.call(
    optimize_thresholds,
    c(list(taxonomy = sorted_tax, refseq = sorted_seq), common_args)
  )
  from_shuffled <- do.call(
    optimize_thresholds,
    c(list(taxonomy = shuffled_tax, refseq = shuffled_seq), common_args)
  )

  # Row identity is (rank, superrank, supertaxon, measure).
  key_cols <- c("rank", "superrank", "supertaxon", "measure")
  from_sorted <- from_sorted[
    do.call(order, from_sorted[key_cols]),
    ,
    drop = FALSE
  ]
  from_shuffled <- from_shuffled[
    do.call(order, from_shuffled[key_cols]),
    ,
    drop = FALSE
  ]
  rownames(from_sorted) <- NULL
  rownames(from_shuffled) <- NULL
  testthat::expect_equal(from_shuffled, from_sorted)
})
