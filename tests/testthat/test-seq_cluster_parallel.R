testthat::test_that("seq_cluster handles one-sequence inputs across methods", {
  threshold_cfg <- threshold_set(c(0.1, 0.2))
  seq_one <- c("ACGTACGTACGT")
  names(seq_one) <- "seq1"

  for (method in c("wfa2", "edlib", "hamming")) {
    for (threads in c(1L, 4L)) {
      msg <- sprintf("method=%s threads=%s", method, threads)
      out <- testthat::expect_no_error(
        seq_cluster(
          seq = seq_one,
          dist_config = dist_config(method = method),
          threshold_config = threshold_cfg,
          clust_config = clust_tree(),
          parallel_config = parallel_concurrent(threads),
          output_type = "matrix",
          which = TRUE,
          verbose = FALSE
        )
      )
      testthat::expect_true(is.matrix(out), info = msg)
      testthat::expect_equal(nrow(out), 2L, info = msg)
      testthat::expect_equal(ncol(out), 1L, info = msg)
    }
  }
})

testthat::test_that("seq_cluster handles zero-sequence inputs across methods", {
  threshold_cfg <- threshold_set(c(0.1, 0.2))

  for (method in c("wfa2", "edlib", "hamming")) {
    for (threads in c(1L, 4L)) {
      msg <- sprintf("method=%s threads=%s", method, threads)
      out <- testthat::expect_no_error(
        seq_cluster(
          seq = character(0),
          seq_id = character(0),
          dist_config = dist_config(method = method),
          threshold_config = threshold_cfg,
          clust_config = clust_tree(),
          parallel_config = parallel_concurrent(threads),
          output_type = "matrix",
          which = TRUE,
          verbose = FALSE
        )
      )
      testthat::expect_true(is.matrix(out), info = msg)
      testthat::expect_equal(nrow(out), 2L, info = msg)
      testthat::expect_equal(ncol(out), 0L, info = msg)
    }
  }
})

testthat::test_that("seq_cluster tiny two-sequence inputs match serial in parallel", {
  seq_two <- c("ACGTACGT", "ACGTTCGT")
  names(seq_two) <- c("seq1", "seq2")
  threshold_cfg <- threshold_set(c(0.1, 0.3))

  serial <- testthat::expect_no_error(
    seq_cluster(
      seq = seq_two,
      dist_config = dist_config(method = "edlib"),
      threshold_config = threshold_cfg,
      clust_config = clust_tree(),
      parallel_config = parallel_concurrent(1L),
      output_type = "matrix",
      verbose = FALSE
    )
  )
  parallel <- testthat::expect_no_error(
    seq_cluster(
      seq = seq_two,
      dist_config = dist_config(method = "edlib"),
      threshold_config = threshold_cfg,
      clust_config = clust_tree(),
      parallel_config = parallel_concurrent(4L),
      output_type = "matrix",
      verbose = FALSE
    )
  )
  testthat::expect_equal(parallel, serial)
})
