testthat::test_that("estimate_subset_memory_mb scales merge by tile shards", {
  n_seq <- 10000L
  n_thresholds <- 40L
  threads <- 2L
  n_tiles <- ceiling(0.5 * (sqrt(1 + 8 * threads) - 1))
  method_scale <- 1 + threads * (2 / n_tiles)
  expected <- max(
    1,
    n_seq * (n_thresholds + 1) * 0.06 * method_scale / 1024
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
  n_tiles_four <- ceiling(0.5 * (sqrt(1 + 8 * 4) - 1))
  expected_four <- max(
    1,
    n_seq * (n_thresholds + 1) * 0.06 * (1 + 4 * (2 / n_tiles_four)) / 1024
  )
  testthat::expect_equal(four, expected_four)
})
