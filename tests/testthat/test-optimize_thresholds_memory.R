testthat::test_that("merge memory scale reflects tile-local shards", {
  n_seq <- 1000L
  n_thresholds <- 401L
  merge8 <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    n_thresholds,
    "tree",
    "merge",
    8L
  )
  concurrent1 <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    n_thresholds,
    "tree",
    "concurrent",
    1L
  )
  # merge(8): parent + 8 * (2/4) tile copies ~= 5x one subset copy
  one_copy <- n_seq * (n_thresholds + 1) * 0.06 / 1024
  testthat::expect_equal(merge8, one_copy * (1 + 8 * (2 / 4)), tolerance = 0.01)
  testthat::expect_lt(merge8, one_copy * 8)
  testthat::expect_gt(merge8, concurrent1)
})
