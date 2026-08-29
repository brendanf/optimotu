# SPDX-CopyrightText: (c) 2025 Brendan Furneaux
# SPDX-License-Identifier: MIT

same_cluster_partition <- function(a, b) {
  cols <- intersect(colnames(a), colnames(b))
  testthat::expect_gt(length(cols), 0L)
  a <- a[, cols, drop = FALSE]
  b <- b[, cols, drop = FALSE]
  testthat::expect_equal(dim(a), dim(b))
  mismatches <- integer()
  for (r in seq_len(nrow(a))) {
    fa <- as.integer(factor(a[r, ], levels = unique(a[r, ])))
    fb <- as.integer(factor(b[r, ], levels = unique(b[r, ])))
    if (!identical(fa, fb)) {
      mismatches <- c(mismatches, r)
    }
  }
  testthat::expect_equal(
    mismatches,
    integer(),
    info = paste(
      "differing rows:",
      paste(rownames(a)[mismatches], collapse = ",")
    )
  )
  invisible(TRUE)
}

make_stream_fixture <- function(n_per = 8L, seq_len = 80L, n_mut = 3L) {
  set.seed(2)
  alphabet <- c("A", "C", "G", "T")
  make_seqs <- function(base, n) {
    vapply(
      seq_len(n),
      function(i) {
        s <- strsplit(base, "", fixed = TRUE)[[1]]
        pos <- sample.int(length(s), n_mut)
        s[pos] <- sample(alphabet, n_mut, replace = TRUE)
        paste(s, collapse = "")
      },
      character(1)
    )
  }
  bases <- replicate(
    2,
    paste(sample(alphabet, seq_len, replace = TRUE), collapse = ""),
    simplify = TRUE
  )
  seq <- unlist(lapply(bases, make_seqs, n = n_per), use.names = FALSE)
  names(seq) <- sprintf("s%03d", seq_along(seq))
  taxonomy <- data.frame(
    seq_id = names(seq),
    kingdom = "k1",
    phylum = rep(c("p1", "p2"), each = n_per),
    class = rep(rep(c("c1", "c2"), each = n_per / 2L), 2),
    stringsAsFactors = FALSE
  )
  list(seq = seq, taxonomy = taxonomy)
}

score_and_pick <- function(k, c, thresh_vals, measures, threads = 1L) {
  k2 <- k
  rownames(k2) <- as.character(seq_len(nrow(k2)))
  df <- calculate_cluster_measures(k2, c, threads, measures)
  confusion <- intersect(measures, c("MCC", "RI", "ARI", "FMI"))
  meas_order <- c(
    confusion,
    intersect(c("MI", "AMI"), measures),
    intersect("FM", measures)
  )
  rows <- lapply(meas_order, function(meas) {
    values <- df$value[as.character(df$measure) == meas]
    best_value <- suppressWarnings(max(values, na.rm = TRUE))
    hits <- which(values == best_value)
    idx <- if (length(hits) == 0L) {
      NA_integer_
    } else {
      strict_median(hits)
    }
    data.frame(
      measure = meas,
      threshold = if (is.na(idx)) NA_real_ else thresh_vals[idx],
      value = best_value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

matrix_optima <- function(
  seq,
  which,
  true_partitions,
  dist_config,
  threshold_config,
  clust_config,
  parallel_config,
  measures
) {
  clust <- seq_cluster(
    seq = seq,
    dist_config = dist_config,
    threshold_config = threshold_config,
    clust_config = clust_config,
    parallel_config = parallel_config,
    output_type = "matrix",
    which = which,
    verbose = FALSE
  )
  thresh_vals <- as.numeric(rownames(clust[[1]]))
  out <- Map(
    function(k, c) {
      score_and_pick(
        k,
        c,
        thresh_vals,
        measures,
        parallel_config$threads
      )
    },
    clust,
    true_partitions
  )
  do.call(rbind, out)
}

testthat::test_that("write_threshold_row matches write_to_matrix partitions", {
  fx <- make_stream_fixture()
  ts <- summarize_by_rank(fx$taxonomy, c("kingdom", "phylum", "class"))
  which <- ts$seq_id
  dist_config <- dist_hamming()
  threshold_config <- threshold_uniform(0, 0.2, 0.02)
  parallel_config <- parallel_concurrent(1L)
  for (clust_config in list(clust_tree(), clust_slink(), clust_matrix())) {
    via_mat <- seq_cluster(
      seq = fx$seq,
      dist_config = dist_config,
      threshold_config = threshold_config,
      clust_config = clust_config,
      parallel_config = parallel_config,
      output_type = "matrix",
      which = which,
      verbose = FALSE
    )
    via_rows <- seq_cluster_multi_via_rows(
      seq = fx$seq,
      which = which,
      dist_config = dist_config,
      threshold_config = threshold_config,
      clust_config = clust_config,
      parallel_config = parallel_config
    )
    testthat::expect_equal(length(via_rows), length(via_mat))
    for (i in seq_along(via_mat)) {
      same_cluster_partition(via_mat[[i]], via_rows[[i]])
    }
  }
})

testthat::test_that("streaming optima match matrix find_best_threshold", {
  fx <- make_stream_fixture()
  ts <- summarize_by_rank(fx$taxonomy, c("kingdom", "phylum", "class"))
  ts <- ts[ts$n_seq >= 4L & ts$n_taxa >= 2L, ]
  measures <- c("MCC", "ARI", "FM")
  dist_config <- dist_hamming()
  threshold_config <- threshold_uniform(0, 0.2, 0.02)
  parallel_config <- parallel_concurrent(1L)
  map <- clustering_threshold_map(threshold_config)
  for (clust_config in list(clust_tree(), clust_slink())) {
    streamed <- seq_cluster_multi_best_threshold(
      seq = fx$seq,
      which = ts$seq_id,
      dist_config = dist_config,
      threshold_config = threshold_config,
      clust_config = clust_config,
      parallel_config = parallel_config,
      true_partitions = ts$true_partition,
      measures = measures,
      thresholds = map$thresholds,
      threshold_order = map$threshold_order
    )
    expected <- matrix_optima(
      seq = fx$seq,
      which = ts$seq_id,
      true_partitions = ts$true_partition,
      dist_config = dist_config,
      threshold_config = threshold_config,
      clust_config = clust_config,
      parallel_config = parallel_config,
      measures = measures
    )
    got <- do.call(rbind, streamed)
    testthat::expect_equal(got$measure, expected$measure)
    testthat::expect_equal(got$threshold, expected$threshold)
    testthat::expect_equal(got$value, expected$value, tolerance = 1e-10)
  }
})

testthat::test_that("optimize_thresholds returns one row per subset and measure", {
  fx <- make_stream_fixture()
  measures <- c("MCC", "ARI")
  out <- optimize_thresholds(
    taxonomy = fx$taxonomy,
    refseq = fx$seq,
    ranks = c("kingdom", "phylum", "class"),
    dist_config = dist_hamming(),
    threshold_config = threshold_uniform(0, 0.2, 0.02),
    clust_config = clust_tree(),
    parallel_config = parallel_concurrent(1L),
    min_taxa = 2L,
    min_refseq = 4L,
    measures = measures,
    verbose = FALSE
  )
  testthat::expect_true(all(
    c("rank", "superrank", "supertaxon", "measure", "threshold", "value") %in%
      names(out)
  ))
  testthat::expect_gt(nrow(out), 0L)
  testthat::expect_true(all(out$measure %in% measures))
})

testthat::test_that("streaming optima handle duplicate thresholds and shuffled which", {
  fx <- make_stream_fixture()
  ts <- summarize_by_rank(fx$taxonomy, c("kingdom", "phylum", "class"))
  ts <- ts[ts$n_seq >= 4L & ts$n_taxa >= 2L, ]
  which <- lapply(ts$seq_id, function(ids) sample(ids))
  true_partitions <- Map(
    function(ids, part, orig) {
      part[match(ids, orig)]
    },
    which,
    ts$true_partition,
    ts$seq_id
  )
  dist_config <- dist_hamming()
  threshold_config <- threshold_set(c(0.02, 0.02, 0.04, 0.08, 0.12))
  clust_config <- clust_tree()
  parallel_config <- parallel_concurrent(1L)
  measures <- c("MCC", "RI")
  map <- clustering_threshold_map(threshold_config)
  streamed <- seq_cluster_multi_best_threshold(
    seq = fx$seq,
    which = which,
    dist_config = dist_config,
    threshold_config = threshold_config,
    clust_config = clust_config,
    parallel_config = parallel_config,
    true_partitions = true_partitions,
    measures = measures,
    thresholds = map$thresholds,
    threshold_order = map$threshold_order
  )
  expected <- matrix_optima(
    seq = fx$seq,
    which = which,
    true_partitions = true_partitions,
    dist_config = dist_config,
    threshold_config = threshold_config,
    clust_config = clust_config,
    parallel_config = parallel_config,
    measures = measures
  )
  got <- do.call(rbind, streamed)
  testthat::expect_equal(got$measure, expected$measure)
  testthat::expect_equal(got$threshold, expected$threshold)
  testthat::expect_equal(got$value, expected$value, tolerance = 1e-10)
})
