testthat::test_that("seq_cluster supports overpartitioned subproblems", {
  seq <- c(
    "ACGTACGTACGT",
    "ACGTACGTTCGT",
    "ACGTTCGTACGT",
    "ACGTTCGTTCGT"
  )
  names(seq) <- paste0("s", seq_along(seq))
  which <- list(
    a = names(seq)[1:3],
    b = names(seq)[2:4]
  )

  baseline <- seq_cluster(
    seq = seq,
    dist_config = dist_edlib(),
    threshold_config = threshold_set(c(0.05, 0.10)),
    clust_config = clust_tree(),
    parallel_config = parallel_merge(threads = 2L),
    output_type = "matrix",
    which = which,
    verbose = FALSE
  )

  overpartitioned_cfg <- parallel_merge(threads = 2L)
  overpartitioned_cfg$subproblems <- 8L
  with_subproblems <- seq_cluster(
    seq = seq,
    dist_config = dist_edlib(),
    threshold_config = threshold_set(c(0.05, 0.10)),
    clust_config = clust_tree(),
    parallel_config = overpartitioned_cfg,
    output_type = "matrix",
    which = which,
    verbose = FALSE
  )

  testthat::expect_equal(with_subproblems, baseline)
})

testthat::test_that("optimize_thresholds budgeted mode preserves output", {
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

  baseline <- suppressWarnings(optimize_thresholds(
    taxonomy = taxonomy,
    refseq = seq,
    ranks = c("kingdom", "phylum", "class"),
    dist_config = dist_edlib(),
    threshold_config = threshold_set(c(0.05, 0.10)),
    clust_config = clust_tree(),
    parallel_config = parallel_merge(threads = 2L),
    min_taxa = 2L,
    min_refseq = 4L,
    measures = c("MCC", "FM"),
    verbose = FALSE
  ))

  budgeted <- suppressWarnings(optimize_thresholds(
    taxonomy = taxonomy,
    refseq = seq,
    ranks = c("kingdom", "phylum", "class"),
    dist_config = dist_edlib(),
    threshold_config = threshold_set(c(0.05, 0.10)),
    clust_config = clust_tree(),
    parallel_config = parallel_merge(threads = 2L),
    min_taxa = 2L,
    min_refseq = 4L,
    measures = c("MCC", "FM"),
    clustering_memory_budget_mb = 64,
    retry_on_memory_exhaustion = TRUE,
    retry_reduce_threads = TRUE,
    retry_overpartition = TRUE,
    verbose = FALSE
  ))

  testthat::expect_equal(budgeted, baseline)
})

testthat::test_that("memory-budget helpers plan, split, and detect errors", {
  testset_select <- data.frame(
    n_seq = c(90L, 40L, 40L, 20L, 20L),
    supertaxon = c("k1", "p1", "p2", "p1", "p2"),
    stringsAsFactors = FALSE
  )
  idx <- seq_len(nrow(testset_select))
  tc <- threshold_set(c(0.05, 0.10, 0.20))
  cc <- clust_tree()
  pc <- parallel_merge(threads = 4L)

  testthat::expect_true(
    optimotu:::is_memory_budget_error(
      simpleError("clustering memory budget exceeded (budget=1.00 MB)")
    )
  )
  testthat::expect_false(
    optimotu:::is_memory_budget_error(simpleError("out of memory"))
  )

  unbudgeted <- optimotu:::make_initial_batches(
    idx,
    testset_select,
    tc,
    cc,
    pc,
    NULL
  )
  testthat::expect_equal(unbudgeted, list(idx))

  tight <- optimotu:::make_initial_batches(
    idx,
    testset_select,
    tc,
    cc,
    pc,
    2
  )
  testthat::expect_gt(length(tight), 1L)
  testthat::expect_equal(sort(unlist(tight)), idx)
  testthat::expect_equal(
    tight,
    optimotu:::make_initial_batches(idx, testset_select, tc, cc, pc, 2)
  )

  overlap <- optimotu:::split_batch_indices(idx, testset_select, "overlap")
  testthat::expect_gt(length(overlap$left), 0L)
  testthat::expect_gt(length(overlap$right), 0L)
  testthat::expect_equal(
    sort(c(overlap$left, overlap$right)),
    idx
  )
  # Overlap grouping keeps the two p1 subsets together.
  testthat::expect_true(
    all(c(2L, 4L) %in% overlap$left) || all(c(2L, 4L) %in% overlap$right)
  )
  testthat::expect_equal(
    overlap,
    optimotu:::split_batch_indices(idx, testset_select, "overlap")
  )

  balanced <- optimotu:::split_batch_indices(idx, testset_select, "balanced")
  testthat::expect_equal(
    sort(c(balanced$left, balanced$right)),
    idx
  )

  reduced <- optimotu:::set_parallel_threads(pc, 2L)
  testthat::expect_equal(reduced$threads, 2L)
  testthat::expect_equal(reduced$method, "merge")
  overpart <- optimotu:::set_parallel_subproblems(pc, 16L)
  testthat::expect_equal(overpart$subproblems, 16L)
})

testthat::test_that("seq_cluster errors when the clustering budget is exceeded", {
  set.seed(1)
  n <- 150L
  seq <- vapply(
    seq_len(n),
    function(i) {
      paste(sample(c("A", "C", "G", "T"), 60L, replace = TRUE), collapse = "")
    },
    character(1)
  )
  names(seq) <- sprintf("s%03d", seq_along(seq))
  which <- lapply(seq_len(12L), function(i) {
    start <- (i - 1L) * 10L + 1L
    names(seq)[seq.int(start, min(n, start + 39L))]
  })
  names(which) <- paste0("w", seq_along(which))

  err <- tryCatch(
    seq_cluster(
      seq = seq,
      dist_config = dist_hamming(),
      threshold_config = threshold_uniform(0, 0.4, 0.01),
      clust_config = clust_tree(),
      parallel_config = parallel_merge(threads = 2L),
      output_type = "matrix",
      which = which,
      clustering_memory_budget_mb = 0.25,
      verbose = FALSE
    ),
    error = identity
  )
  testthat::expect_s3_class(err, "error")
  testthat::expect_true(optimotu:::is_memory_budget_error(err))
})

testthat::test_that("tight-budget optimize_thresholds matches the unbudgeted result", {
  set.seed(1)
  n_per <- 20L
  alphabet <- c("A", "C", "G", "T")
  make_seqs <- function(base, n) {
    vapply(
      seq_len(n),
      function(i) {
        s <- strsplit(base, "", fixed = TRUE)[[1]]
        pos <- sample.int(length(s), 3L)
        s[pos] <- sample(alphabet, 3L, replace = TRUE)
        paste(s, collapse = "")
      },
      character(1)
    )
  }
  bases <- replicate(
    3,
    paste(sample(alphabet, 80L, replace = TRUE), collapse = ""),
    simplify = TRUE
  )
  seq <- unlist(lapply(bases, make_seqs, n = n_per), use.names = FALSE)
  names(seq) <- sprintf("s%03d", seq_along(seq))
  taxonomy <- data.frame(
    seq_id = names(seq),
    kingdom = "k1",
    phylum = rep(c("p1", "p2", "p3"), each = n_per),
    class = rep(rep(c("c1", "c2"), each = n_per / 2L), 3),
    stringsAsFactors = FALSE
  )

  args <- list(
    taxonomy = taxonomy,
    refseq = seq,
    ranks = c("kingdom", "phylum", "class"),
    dist_config = dist_hamming(),
    threshold_config = threshold_uniform(0, 0.4, 0.01),
    clust_config = clust_tree(),
    parallel_config = parallel_merge(threads = 2L),
    min_taxa = 2L,
    min_refseq = 4L,
    measures = c("MCC", "FM"),
    verbose = FALSE
  )
  baseline <- do.call(optimize_thresholds, args)
  budgeted <- do.call(
    optimize_thresholds,
    c(args, list(clustering_memory_budget_mb = 1))
  )
  testthat::expect_equal(budgeted, baseline)
})
