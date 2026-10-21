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

# Square all-pairs side length. AllPairGenerator builds n_blocks^2 tiles.
merge_n_blocks <- function(threads) {
  ceiling(sqrt(threads))
}

merge_n_slots <- function(parallel_method = "merge", threads = 2L) {
  if (parallel_method == "merge" && threads > 1) {
    n_blocks <- merge_n_blocks(threads)
    return(as.integer(n_blocks * n_blocks))
  }
  if (parallel_method == "hierarchical") {
    return(as.integer(max(threads, 1L)))
  }
  as.integer(max(threads, 1L))
}

kruskal_parent_bytes <- function(
  n_seq,
  parallel_method = "merge",
  threads = 2L
) {
  n_slots <- merge_n_slots(parallel_method, threads)
  n_pad <- 1
  while (n_pad < n_slots) {
    n_pad <- n_pad * 2
  }
  conc <- max(1L, n_slots %/% 2L)
  2 *
    n_pad *
    32 +
    n_slots * n_seq * 16 +
    conc * (n_seq * 8 + 3 * n_seq * 16) +
    n_seq * 16
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
      kruskal_parent_bytes(n_seq, parallel_method, threads)
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
  n_blocks <- merge_n_blocks(threads)
  1 + threads * (2 / n_blocks)
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
  # 2 threads -> n_blocks=2, n_slots=4, n_pad=4, conc=2.
  testthat::expect_equal(
    kruskal_parent_bytes(n_seq, "merge", threads),
    2 * 4 * 32 + 4 * n_seq * 16 + 2 * (n_seq * 8 + 3 * n_seq * 16) + n_seq * 16
  )
})

testthat::test_that("merge slink parent uses Kruskal tournament formula", {
  n_seq <- 10000L
  n_thresholds <- 40L
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
  # 120 threads -> n_blocks=11, n_slots=121.
  testthat::expect_equal(
    kruskal_parent_bytes(n_seq, "merge", 120L),
    {
      n_slots <- 11L * 11L
      n_pad <- 128L
      conc <- n_slots %/% 2L
      2 *
        n_pad *
        32 +
        n_slots * n_seq * 16 +
        conc * (n_seq * 8 + 3 * n_seq * 16) +
        n_seq * 16
    }
  )
  testthat::expect_gt(
    kruskal_parent_bytes(n_seq, "merge", 120L),
    n_seq * 8 * 32
  )
})

testthat::test_that("estimate_subset_memory_mb merge scale grows with threads", {
  n_seq <- 10000L
  n_thresholds <- 40L
  # 2 threads: n_blocks=2, scale=3. 4 threads: n_blocks=2, scale=5.
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

testthat::test_that("incremental batch planner matches unique(unlist) unions", {
  # Overlapping character IDs and explicit seq_index must agree, and batch
  # cuts must match a naive unique(unlist()) planner.
  testset <- data.frame(
    n_seq = c(10L, 12L, 8L, 15L),
    stringsAsFactors = FALSE
  )
  testset$seq_id <- list(
    paste0("s", 1:10),
    paste0("s", 5:16),
    paste0("s", 14:21),
    paste0("s", 1:15)
  )
  testset$seq_index <- list(
    1:10,
    5:16,
    14:21,
    1:15
  )
  threshold_config <- threshold_uniform(0, 0.1, 0.05)
  clust_config <- clust_slink()
  parallel_config <- parallel_concurrent(1L)
  ordered_idx <- seq_len(nrow(testset))

  naive_batches <- function(budget_mb) {
    batches <- list()
    current <- integer(0)
    for (idx in ordered_idx) {
      candidate <- c(current, idx)
      est <- optimotu:::estimate_batch_memory_mb(
        candidate,
        testset,
        threshold_config,
        clust_config,
        parallel_config,
        dist_method = "hamming",
        mean_seq_len = 400
      )
      if (length(current) > 0L && est > budget_mb) {
        batches[[length(batches) + 1L]] <- current
        current <- idx
      } else {
        current <- candidate
      }
    }
    batches[[length(batches) + 1L]] <- current
    batches
  }

  for (budget in c(0.05, 0.2, 1, 100)) {
    got <- optimotu:::make_initial_batches(
      ordered_idx,
      testset,
      threshold_config,
      clust_config,
      parallel_config,
      budget,
      dist_method = "hamming",
      mean_seq_len = 400
    )
    testthat::expect_equal(got, naive_batches(budget), info = budget)
  }

  # seq_id-only fallback agrees with seq_index path.
  testset_ids_only <- testset
  testset_ids_only$seq_index <- NULL
  testthat::expect_equal(
    optimotu:::make_initial_batches(
      ordered_idx,
      testset,
      threshold_config,
      clust_config,
      parallel_config,
      0.2,
      dist_method = "hamming",
      mean_seq_len = 400
    ),
    optimotu:::make_initial_batches(
      ordered_idx,
      testset_ids_only,
      threshold_config,
      clust_config,
      parallel_config,
      0.2,
      dist_method = "hamming",
      mean_seq_len = 400
    )
  )
})

testthat::test_that("estimate_batch_memory_mb adds Hamming packed union bytes", {
  testset <- data.frame(
    n_seq = c(10L, 12L, 8L),
    stringsAsFactors = FALSE
  )
  testset$seq_id <- list(
    paste0("s", 1:10),
    paste0("s", 5:16),
    paste0("s", 14:21)
  )
  threshold_config <- threshold_uniform(0, 0.1, 0.05)
  clust_config <- clust_slink()
  parallel_config <- parallel_concurrent(1L)
  without <- optimotu:::estimate_batch_memory_mb(
    batch_idx = 1:2,
    testset_select = testset,
    threshold_config = threshold_config,
    clust_config = clust_config,
    parallel_config = parallel_config
  )
  with_ham <- optimotu:::estimate_batch_memory_mb(
    batch_idx = 1:2,
    testset_select = testset,
    threshold_config = threshold_config,
    clust_config = clust_config,
    parallel_config = parallel_config,
    dist_method = "hamming",
    mean_seq_len = 400
  )
  n_union <- length(unique(unlist(testset$seq_id[1:2], use.names = FALSE)))
  packed <- optimotu:::estimate_hamming_packed_mb(n_union, 400)
  # Non-Hamming path already includes SequenceView allowance.
  testthat::expect_equal(with_ham, without + packed)
  testthat::expect_gt(with_ham, without)
})
