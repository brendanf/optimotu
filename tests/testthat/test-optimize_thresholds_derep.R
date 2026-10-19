# SPDX-CopyrightText: (c) 2026 Brendan Furneaux
# SPDX-License-Identifier: MIT

make_dup_fixture <- function() {
  set.seed(11)
  alphabet <- c("A", "C", "G", "T")
  bases <- replicate(
    4,
    paste(sample(alphabet, 60L, replace = TRUE), collapse = ""),
    simplify = TRUE
  )
  seq <- c(
    bases[1],
    bases[1],
    bases[1],
    bases[2],
    bases[3],
    bases[3],
    bases[4]
  )
  names(seq) <- sprintf("s%02d", seq_along(seq))
  taxonomy <- data.frame(
    seq_id = names(seq),
    kingdom = "k1",
    phylum = c("p1", "p1", "p1", "p1", "p2", "p2", "p2"),
    class = c("c1", "c1", "c1", "c2", "c3", "c3", "c4"),
    stringsAsFactors = FALSE
  )
  list(seq = seq, taxonomy = taxonomy)
}

make_conflict_fixture <- function() {
  fx <- make_dup_fixture()
  # Same bytes as s01, different class: must not collapse with c1 copies.
  fx$seq <- c(fx$seq, fx$seq[["s01"]])
  names(fx$seq)[length(fx$seq)] <- "s08"
  fx$taxonomy <- rbind(
    fx$taxonomy,
    data.frame(
      seq_id = "s08",
      kingdom = "k1",
      phylum = "p1",
      class = "c2",
      stringsAsFactors = FALSE
    )
  )
  fx
}

oracle_optima <- function(
  taxonomy,
  refseq,
  ranks,
  dist_config,
  threshold_config,
  clust_config,
  parallel_config,
  min_taxa,
  min_refseq,
  measures,
  ...
) {
  seq <- refseq
  taxonomy <- taxonomy[match(names(seq), taxonomy$seq_id), , drop = FALSE]
  ts <- summarize_by_rank(taxonomy, ranks, "seq_id")
  ts <- ts[
    (ts$n_seq >= min_refseq & ts$n_taxa >= min_taxa) |
      ts$supertaxon %in% taxonomy[[ranks[1]]],
  ]
  clust <- seq_cluster(
    seq = seq,
    dist_config = dist_config,
    threshold_config = threshold_config,
    clust_config = clust_config,
    parallel_config = parallel_config,
    output_type = "matrix",
    which = ts$seq_id,
    verbose = FALSE
  )
  rows <- Map(
    function(k, c, i) {
      best <- suppressWarnings(find_best_threshold(
        k,
        c,
        threads = parallel_config$threads,
        measures = measures
      ))
      cbind(
        rank = ts$rank[i],
        superrank = ts$superrank[i],
        supertaxon = ts$supertaxon[i],
        best,
        stringsAsFactors = FALSE
      )
    },
    clust,
    ts$true_partition,
    seq_along(clust)
  )
  do.call(rbind, rows)
}

sort_optima <- function(df) {
  df[order(df$rank, df$superrank, df$supertaxon, df$measure), , drop = FALSE]
}

expect_optima_equal <- function(got, expected) {
  got <- sort_optima(got)
  expected <- sort_optima(expected)
  testthat::expect_equal(
    got[, c("rank", "superrank", "supertaxon", "measure")],
    expected[, c("rank", "superrank", "supertaxon", "measure")],
    ignore_attr = TRUE
  )
  testthat::expect_equal(got$threshold, expected$threshold)
  testthat::expect_equal(got$value, expected$value, tolerance = 1e-10)
}

common_args <- function(
  fx,
  threshold_config = threshold_uniform(0, 0.2, 0.02),
  measures = c("MCC", "ARI", "FM", "AMI")
) {
  list(
    taxonomy = fx$taxonomy,
    refseq = fx$seq,
    ranks = c("kingdom", "phylum", "class"),
    dist_config = dist_hamming(),
    threshold_config = threshold_config,
    clust_config = clust_tree(),
    parallel_config = parallel_concurrent(threads = 1L),
    min_taxa = 2L,
    min_refseq = 3L,
    measures = measures,
    verbose = FALSE
  )
}

testthat::test_that("group_identical_seq_tax keeps first-seen representatives", {
  fx <- make_dup_fixture()
  grouped <- group_identical_seq_tax(
    fx$seq,
    fx$taxonomy[c("kingdom", "phylum", "class")]
  )
  testthat::expect_equal(grouped$representative, c(1L, 4L, 5L, 7L))
  testthat::expect_equal(grouped$map, c(1L, 1L, 1L, 2L, 3L, 3L, 4L))
})

testthat::test_that("identical sequence with different taxonomy is not grouped", {
  fx <- make_conflict_fixture()
  grouped <- group_identical_seq_tax(
    fx$seq,
    fx$taxonomy[c("kingdom", "phylum", "class")]
  )
  testthat::expect_equal(length(grouped$representative), 5L)
  testthat::expect_true(8L %in% grouped$representative)
  testthat::expect_false(grouped$map[8] == grouped$map[1])
})

testthat::test_that("optimize_thresholds matches unreduced oracle with duplicates", {
  fx <- make_dup_fixture()
  args <- common_args(fx)
  got <- do.call(optimize_thresholds, args)
  expected <- do.call(oracle_optima, args)
  expect_optima_equal(got, expected)
})

testthat::test_that("unweighted representatives differ from weighted oracle", {
  fx <- make_dup_fixture()
  args <- common_args(fx, measures = c("MCC", "ARI", "FM"))
  taxonomy <- fx$taxonomy[match(names(fx$seq), fx$taxonomy$seq_id), ]
  ts <- summarize_by_rank(taxonomy, args$ranks, "seq_id")
  ts <- ts[
    (ts$n_seq >= args$min_refseq & ts$n_taxa >= args$min_taxa) |
      ts$supertaxon %in% taxonomy[[args$ranks[1]]],
  ]
  grouped <- group_identical_seq_tax(fx$seq, taxonomy[args$ranks])
  seq_rep <- fx$seq[grouped$representative]
  remapped <- remap_testset_representatives(
    ts$seq_index,
    ts$true_partition,
    grouped$map,
    names(seq_rep)
  )
  map <- clustering_threshold_map(args$threshold_config)
  which_idx <- lapply(remapped$seq_index, function(i) as.integer(i - 1L))
  unweighted <- seq_cluster_multi_best_threshold(
    seq = seq_rep,
    which = which_idx,
    dist_config = args$dist_config,
    threshold_config = args$threshold_config,
    clust_config = args$clust_config,
    parallel_config = args$parallel_config,
    true_partitions = remapped$true_partition,
    measures = args$measures,
    thresholds = map$thresholds,
    threshold_order = map$threshold_order
  )
  unweighted <- do.call(rbind, unweighted)
  expected <- do.call(oracle_optima, args)
  testthat::expect_false(isTRUE(all.equal(
    unweighted$value,
    expected$value,
    tolerance = 1e-10
  )))
})

testthat::test_that("optimize_thresholds is unchanged with no duplicates", {
  fx <- make_dup_fixture()
  keep <- !duplicated(fx$seq)
  fx$seq <- fx$seq[keep]
  fx$taxonomy <- fx$taxonomy[keep, , drop = FALSE]
  args <- common_args(fx)
  got <- do.call(optimize_thresholds, args)
  expected <- do.call(oracle_optima, args)
  expect_optima_equal(got, expected)
})

testthat::test_that("conflicting taxonomy copies are not collapsed", {
  fx <- make_conflict_fixture()
  args <- common_args(fx)
  got <- do.call(optimize_thresholds, args)
  expected <- do.call(oracle_optima, args)
  expect_optima_equal(got, expected)
})

testthat::test_that("threshold grid of only 0 does not collapse copies", {
  fx <- make_dup_fixture()
  args <- common_args(fx, threshold_config = threshold_set(0))
  testthat::expect_false(
    thresholds_allow_identical_collapse(args$threshold_config)
  )
  got <- do.call(optimize_thresholds, args)
  expected <- do.call(oracle_optima, args)
  expect_optima_equal(got, expected)
})
