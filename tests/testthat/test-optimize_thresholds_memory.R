algo_bytes <- function(n_seq, n_thresholds, clust_method) {
  int_bytes <- 4
  u64_bytes <- 8
  switch(
    clust_method,
    tree = n_seq * u64_bytes * 32,
    matrix = n_seq * n_thresholds * int_bytes + n_seq * int_bytes * 8,
    index = n_seq * n_thresholds * int_bytes + n_seq * int_bytes * 12,
    n_seq * n_thresholds * int_bytes
  )
}

slink_leaf_bytes <- function(n_seq) {
  n_seq * 8 * 4
}

slink_parent_bytes <- function(n_seq, parallel_method = "merge", threads = 2L) {
  leaf <- slink_leaf_bytes(n_seq)
  if (parallel_method == "merge" && threads > 1) {
    n_tiles <- ceiling(0.5 * (sqrt(1 + 8 * threads) - 1))
    return(leaf + n_seq * n_tiles * 16)
  }
  if (parallel_method == "hierarchical") {
    return(leaf + n_seq * threads * 16)
  }
  leaf + n_seq * max(threads, 1L) * 16
}

estimate_subset_bytes <- function(
  n_seq,
  n_thresholds,
  clust_method,
  parallel_method,
  threads
) {
  method_scale <- merge_scale(threads)
  if (clust_method == "slink") {
    root <- if (
      parallel_method == "concurrent" ||
        (parallel_method == "merge" && threads <= 1)
    ) {
      slink_leaf_bytes(n_seq)
    } else {
      slink_parent_bytes(n_seq, parallel_method, threads)
    }
    child <- slink_leaf_bytes(n_seq)
    scale <- switch(
      parallel_method,
      merge = if (threads <= 1) 1.1 else method_scale,
      concurrent = 1.1,
      hierarchical = max(1, threads / 2),
      1
    )
    return(root + (scale - 1) * child)
  }
  scale <- switch(
    parallel_method,
    merge = if (threads <= 1) 1.1 else method_scale,
    concurrent = 1.1,
    hierarchical = max(1, threads / 2),
    1
  )
  algo_bytes(n_seq, n_thresholds, clust_method) * scale
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
      estimate_subset_bytes(
        n_seq,
        n_thresholds,
        method,
        "concurrent",
        1L
      ) /
        (1024 * 1024)
    )
    testthat::expect_equal(got, expected, info = method)
  }
})

testthat::test_that("one-thread merge matches concurrent estimates", {
  n_seq <- 10000L
  n_thresholds <- 40L
  for (method in c("tree", "slink", "matrix", "index")) {
    merge_one <- optimotu:::estimate_subset_memory_mb(
      n_seq,
      n_thresholds,
      method,
      "merge",
      1L
    )
    concurrent_one <- optimotu:::estimate_subset_memory_mb(
      n_seq,
      n_thresholds,
      method,
      "concurrent",
      1L
    )
    testthat::expect_equal(merge_one, concurrent_one, info = method)
  }
})

testthat::test_that("concurrent slink estimate is smaller than concurrent tree", {
  n_seq <- 10000L
  n_thresholds <- 40L
  tree_mb <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    n_thresholds,
    "tree",
    "concurrent",
    1L
  )
  slink_mb <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    n_thresholds,
    "slink",
    "concurrent",
    1L
  )
  testthat::expect_lt(slink_mb, tree_mb)
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
  result_bytes <- n_seq * n_thresholds * 4
  for (method in c("tree", "slink")) {
    base_bytes <- estimate_subset_bytes(
      n_seq,
      n_thresholds,
      method,
      "concurrent",
      1L
    )
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
    expected <- max(1, (base_bytes + result_bytes) / (1024 * 1024))
    testthat::expect_equal(with_result, expected, info = method)
    testthat::expect_equal(
      without,
      max(1, base_bytes / (1024 * 1024)),
      info = method
    )
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

testthat::test_that("merge slink uses parent root and leaf tile shards", {
  n_seq <- 10000L
  n_thresholds <- 40L
  threads <- 2L
  expected <- max(
    1,
    estimate_subset_bytes(
      n_seq,
      n_thresholds,
      "slink",
      "merge",
      threads
    ) /
      (1024 * 1024)
  )
  got <- optimotu:::estimate_subset_memory_mb(
    n_seq = n_seq,
    n_thresholds = n_thresholds,
    clust_method = "slink",
    parallel_method = "merge",
    threads = threads
  )
  testthat::expect_equal(got, expected)
  tree_merge <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    n_thresholds,
    "tree",
    "merge",
    threads
  )
  testthat::expect_lt(got, tree_merge)
  # Parent cache peak at T=2 is 64n, below tree 256n.
  testthat::expect_equal(
    slink_parent_bytes(n_seq, "merge", threads),
    slink_leaf_bytes(n_seq) + n_seq * 2 * 16
  )
})

testthat::test_that("merge slink parent uses uncapped cache formula", {
  n_seq <- 10000L
  n_thresholds <- 40L
  # T=8 at 32 threads: cache peak 160n < 256n tree.
  slink_32 <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    n_thresholds,
    "slink",
    "merge",
    32L
  )
  tree_32 <- optimotu:::estimate_subset_memory_mb(
    n_seq,
    n_thresholds,
    "tree",
    "merge",
    32L
  )
  testthat::expect_lt(slink_32, tree_32)
  # Above former tree-sized crossover (T > 14), estimate keeps growing.
  n_tiles_120 <- ceiling(0.5 * (sqrt(1 + 8 * 120) - 1))
  testthat::expect_equal(
    slink_parent_bytes(n_seq, "merge", 120L),
    slink_leaf_bytes(n_seq) + n_seq * n_tiles_120 * 16
  )
  testthat::expect_gt(
    slink_parent_bytes(n_seq, "merge", 120L),
    n_seq * 8 * 32
  )
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
