algo_bytes <- function(n_seq, n_thresholds, clust_method) {
  int_bytes <- 4
  u64_bytes <- 8
  switch(
    clust_method,
    tree = n_seq * u64_bytes * 32,
    slink = n_seq * u64_bytes * 24,
    matrix = n_seq * n_thresholds * int_bytes + n_seq * int_bytes * 8,
    index = n_seq * n_thresholds * int_bytes + n_seq * int_bytes * 12,
    n_seq * n_thresholds * int_bytes
  )
}

merge_scale <- function(threads) {
  n_tiles <- ceiling(0.5 * (sqrt(1 + 8 * threads) - 1))
  1 + threads * (2 / n_tiles)
}

testthat::test_that("estimate_subset_memory_mb matches estimate_bytes logic", {
  n_seq <- 10000L
  n_thresholds <- 40L
  for (method in c("tree", "slink", "matrix", "index")) {
    got <- optimotu:::estimate_subset_memory_mb(
      n_seq = n_seq,
      n_thresholds = n_thresholds,
      clust_method = method,
      parallel_method = "concurrent",
      threads = 1L
    )
    expected <- max(
      1,
      algo_bytes(n_seq, n_thresholds, method) * 1.1 / (1024 * 1024)
    )
    testthat::expect_equal(got, expected, info = method)
  }
})

testthat::test_that("tree and slink estimates do not grow with n_thresholds", {
  n_seq <- 10000L
  tree_few <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    10L,
    "tree",
    "concurrent",
    1L
  )
  tree_many <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    400L,
    "tree",
    "concurrent",
    1L
  )
  slink_few <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    10L,
    "slink",
    "concurrent",
    1L
  )
  slink_many <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    400L,
    "slink",
    "concurrent",
    1L
  )
  testthat::expect_equal(tree_few, tree_many)
  testthat::expect_equal(slink_few, slink_many)
  matrix_few <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    10L,
    "matrix",
    "concurrent",
    1L
  )
  matrix_many <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    400L,
    "matrix",
    "concurrent",
    1L
  )
  testthat::expect_gt(matrix_many, matrix_few)
})

testthat::test_that("include_result adds 4nm only for tree and slink", {
  n_seq <- 10000L
  n_thresholds <- 40L
  result_mb <- n_seq * n_thresholds * 4 / (1024 * 1024)
  for (method in c("tree", "slink")) {
    without <- optimotu:::estimate_subset_memory_mb(
      n_seq,
      n_thresholds,
      method,
      "concurrent",
      1L,
      include_result = FALSE
    )
    with_result <- optimotu:::estimate_subset_memory_mb(
      n_seq,
      n_thresholds,
      method,
      "concurrent",
      1L,
      include_result = TRUE
    )
    testthat::expect_equal(with_result, without + result_mb, info = method)
  }
  for (method in c("matrix", "index")) {
    without <- optimotu:::estimate_subset_memory_mb(
      n_seq,
      n_thresholds,
      method,
      "concurrent",
      1L,
      include_result = FALSE
    )
    with_result <- optimotu:::estimate_subset_memory_mb(
      n_seq,
      n_thresholds,
      method,
      "concurrent",
      1L,
      include_result = TRUE
    )
    testthat::expect_equal(with_result, without, info = method)
  }
})

testthat::test_that("estimate_subset_memory_mb scales merge by tile shards", {
  n_seq <- 10000L
  n_thresholds <- 40L
  threads <- 2L
  expected <- max(
    1,
    algo_bytes(n_seq, n_thresholds, "tree") *
      merge_scale(threads) /
      (1024 * 1024)
  )
  got <- optimotu:::estimate_subset_memory_mb(
    n_seq = n_seq,
    n_thresholds = n_thresholds,
    clust_method = "tree",
    parallel_method = "merge",
    threads = threads
  )
  testthat::expect_equal(got, expected)
})

testthat::test_that("estimate_subset_memory_mb merge scale grows with threads", {
  n_seq <- 10000L
  n_thresholds <- 40L
  # threads=1 and threads=2 both use n_tiles such that the scale is 3.
  two <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    n_thresholds,
    "tree",
    "merge",
    2L
  )
  four <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    n_thresholds,
    "tree",
    "merge",
    4L
  )
  testthat::expect_gt(four, two)
  expected_four <- max(
    1,
    algo_bytes(n_seq, n_thresholds, "tree") *
      merge_scale(4L) /
      (1024 * 1024)
  )
  testthat::expect_equal(four, expected_four)
})
