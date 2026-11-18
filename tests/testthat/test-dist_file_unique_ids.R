# SPDX-CopyrightText: (c) 2026 Brendan Furneaux
# SPDX-License-Identifier: MIT

make_unique_id_fixture <- function() {
  set.seed(21)
  alphabet <- c("A", "C", "G", "T")
  seq <- replicate(
    6L,
    paste(sample(alphabet, 40L, replace = TRUE), collapse = ""),
    simplify = TRUE
  )
  names(seq) <- sprintf("s%02d", seq_along(seq))
  taxonomy <- data.frame(
    seq_id = names(seq),
    kingdom = "k1",
    phylum = rep(c("p1", "p2"), each = 3L),
    class = rep(c("c1", "c2", "c3"), 2L),
    stringsAsFactors = FALSE
  )
  list(seq = seq, taxonomy = taxonomy)
}

write_named_dist <- function(seq, path) {
  n <- length(seq)
  rows <- list()
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      d <- mean(
        strsplit(seq[[i]], "", fixed = TRUE)[[1]] !=
          strsplit(seq[[j]], "", fixed = TRUE)[[1]]
      )
      rows[[length(rows) + 1L]] <- data.frame(
        id1 = names(seq)[i],
        id2 = names(seq)[j],
        dist = d,
        stringsAsFactors = FALSE
      )
    }
  }
  write.table(
    do.call(rbind, rows),
    file = path,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  invisible(path)
}

testthat::test_that("optimize_thresholds rejects duplicate taxonomy IDs", {
  fx <- make_unique_id_fixture()
  dup_tax <- rbind(fx$taxonomy, fx$taxonomy[1, , drop = FALSE])
  dup_seq <- c(fx$seq, fx$seq[1])
  names(dup_seq)[length(dup_seq)] <- names(fx$seq)[1]

  testthat::expect_error(
    optimize_thresholds(
      taxonomy = dup_tax,
      refseq = dup_seq,
      ranks = c("kingdom", "phylum", "class"),
      dist_config = dist_hamming(),
      threshold_config = threshold_set(c(0.05, 0.10, 0.20)),
      clust_config = clust_tree(),
      parallel_config = parallel_concurrent(1L),
      min_taxa = 2L,
      min_refseq = 2L,
      measures = "MCC",
      verbose = FALSE
    ),
    "taxonomy\\[\\['seq_id'\\]\\].*Contains duplicated values"
  )

  testthat::expect_no_error(
    optimize_thresholds(
      taxonomy = fx$taxonomy,
      refseq = fx$seq,
      ranks = c("kingdom", "phylum", "class"),
      dist_config = dist_hamming(),
      threshold_config = threshold_set(c(0.05, 0.10, 0.20)),
      clust_config = clust_tree(),
      parallel_config = parallel_concurrent(1L),
      min_taxa = 2L,
      min_refseq = 2L,
      measures = "MCC",
      verbose = FALSE
    )
  )
})

testthat::test_that("seq_cluster_usearch rejects duplicate sequence names before USEARCH", {
  seq <- c(
    s01 = "ACGTACGTACGTACGTACGT",
    s02 = "ACGTACGTACGTACGTACTT",
    s03 = "ACGTACGTACGTACGTTCGT"
  )
  testthat::expect_error(
    seq_cluster_usearch(
      seq = seq,
      seq_id = c("s01", "s01", "s03"),
      threshold_config = threshold_set(c(0.05, 0.10)),
      clust_config = clust_tree(),
      parallel_config = parallel_concurrent(1L),
      usearch = "/nonexistent/usearch"
    ),
    "dist_usearch\\(\\) sequence names.*Contains duplicated values"
  )
})

testthat::test_that("optimize_thresholds rejects duplicate USEARCH names before system2", {
  taxonomy <- data.frame(
    seq_id = c("a", "b"),
    kingdom = "k1",
    phylum = c("p1", "p1"),
    class = c("c1", "c2"),
    stringsAsFactors = FALSE
  )
  # Two FASTA headers share the taxonomy set {a,b}; match() expands taxonomy.
  seq <- c(
    a = "ACGTACGTACGTACGTACGT",
    a = "ACGTACGTACGTACGTACTT",
    b = "ACGTACGTACGTACGTTCGT"
  )
  dist_cfg <- structure(
    list(
      method = "usearch",
      usearch = "/nonexistent/usearch",
      usearch_ncpu = 1L,
      call = quote(dist_usearch())
    ),
    class = "optimotu_dist_config"
  )
  testthat::expect_error(
    optimize_thresholds(
      taxonomy = taxonomy,
      refseq = seq,
      ranks = c("kingdom", "phylum", "class"),
      dist_config = dist_cfg,
      threshold_config = threshold_set(c(0.05, 0.10)),
      clust_config = clust_tree(),
      parallel_config = parallel_concurrent(1L),
      min_taxa = 2L,
      min_refseq = 2L,
      measures = "MCC",
      verbose = FALSE
    ),
    "dist_usearch\\(\\) sequence names.*Contains duplicated values"
  )
})

testthat::test_that("dist_file path rejects duplicate names and accepts unique names", {
  fx <- make_unique_id_fixture()
  dist_path <- tempfile(fileext = ".dist")
  on.exit(unlink(dist_path), add = TRUE)
  write_named_dist(fx$seq, dist_path)

  dup_tax <- rbind(fx$taxonomy, fx$taxonomy[1, , drop = FALSE])
  dup_seq <- c(fx$seq, fx$seq[1])
  names(dup_seq)[length(dup_seq)] <- names(fx$seq)[1]

  testthat::expect_error(
    optimize_thresholds(
      taxonomy = dup_tax,
      refseq = dup_seq,
      ranks = c("kingdom", "phylum", "class"),
      dist_config = dist_file(dist_path, by_name = TRUE),
      threshold_config = threshold_set(c(0.05, 0.10, 0.20)),
      clust_config = clust_tree(),
      parallel_config = parallel_concurrent(1L),
      min_taxa = 2L,
      min_refseq = 2L,
      measures = "MCC",
      verbose = FALSE
    ),
    "taxonomy\\[\\['seq_id'\\]\\].*Contains duplicated values"
  )

  testthat::expect_no_error(
    optimize_thresholds(
      taxonomy = fx$taxonomy,
      refseq = fx$seq,
      ranks = c("kingdom", "phylum", "class"),
      dist_config = dist_file(dist_path, by_name = TRUE),
      threshold_config = threshold_set(c(0.05, 0.10, 0.20)),
      clust_config = clust_tree(),
      parallel_config = parallel_concurrent(1L),
      min_taxa = 2L,
      min_refseq = 2L,
      measures = "MCC",
      verbose = FALSE
    )
  )
})

testthat::test_that("distmx_cluster rejects duplicate sequence names", {
  dist_path <- tempfile(fileext = ".dist")
  on.exit(unlink(dist_path), add = TRUE)
  writeLines(
    c(
      "s01\ts02\t0.1",
      "s01\ts03\t0.2",
      "s02\ts03\t0.15"
    ),
    dist_path
  )
  testthat::expect_error(
    distmx_cluster(
      distmx = dist_path,
      names = c("s01", "s01", "s02"),
      threshold_config = threshold_set(c(0.05, 0.10, 0.20)),
      clust_config = clust_tree(),
      parallel_config = parallel_concurrent(1L),
      by_name = TRUE
    ),
    "duplicate sequence name 's01'"
  )

  testthat::expect_error(
    distmx_cluster(
      distmx = dist_path,
      names = c("s01", "s02", "s03"),
      threshold_config = threshold_set(c(0.05, 0.10, 0.20)),
      clust_config = clust_tree(),
      parallel_config = parallel_concurrent(1L),
      which = list(c("s01", "s01", "s02")),
      by_name = TRUE
    ),
    "sequence 's01' is listed more than once in a subset"
  )
})

testthat::test_that("seq_cluster rejects a repeated member in one which subset", {
  seq <- c(
    s01 = "ACGTACGTACGTACGTACGT",
    s02 = "ACGTACGTACGTACGTACTT",
    s03 = "ACGTACGTACGTACGTTCGT",
    s04 = "ACGTACGTACGTACTTCCGT"
  )
  testthat::expect_error(
    seq_cluster(
      seq = seq,
      dist_config = dist_hamming(),
      threshold_config = threshold_set(c(0.05, 0.10, 0.20)),
      clust_config = clust_tree(),
      parallel_config = parallel_concurrent(1L),
      which = list(c("s01", "s01", "s02", "s03"))
    ),
    "sequence 's01' is listed more than once in a subset"
  )

  # Integer which (0-based) hits sort_subset_indices() rather than the
  # name-based router.
  testthat::expect_error(
    seq_cluster_multi_best_threshold(
      seq = unname(seq),
      which = list(c(0L, 0L, 1L, 2L)),
      dist_config = dist_hamming(),
      threshold_config = threshold_set(c(0.05, 0.10, 0.20)),
      clust_config = clust_tree(),
      parallel_config = parallel_concurrent(1L),
      true_partitions = list(c(1L, 1L, 2L, 3L)),
      measures = "MCC",
      thresholds = c(0.05, 0.10, 0.20),
      threshold_order = 1:3
    ),
    "subset index 0 is listed more than once in a subset"
  )
})
