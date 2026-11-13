test_that("score mode matches identity on substitution-only pairs", {
  a <- "ACGTACGT"
  b <- "ACGTACCT"
  for (method in c("wfa2", "edlib", "ksw2")) {
    for (span in c("global", "extend")) {
      d_id <- align(a, b, method = method, span = span)
      d_sc <- align(
        a,
        b,
        method = method,
        span = span,
        dissimilarity = "score"
      )
      if (method == "ksw2" && span == "extend") {
        # KSW2 identity-extend can stop before a trailing substitution;
        # score mode still counts it.
        expect_equal(d_sc, 0.125)
      } else {
        expect_equal(d_sc, d_id)
      }
    }
  }
  # KSW2 match = 0 maps to mu = 1
  expect_equal(
    align("ACGT", "ACGA", method = "ksw2", match = 0, dissimilarity = "score"),
    align("ACGT", "ACGA", method = "ksw2", match = 0)
  )
})

test_that("extend scores a perfect prefix as 0; global does not", {
  long <- "AAACGCGTCG"
  short <- "AAACGCGT"
  for (method in c("wfa2", "edlib", "ksw2")) {
    expect_equal(
      align(
        long,
        short,
        method = method,
        span = "extend",
        dissimilarity = "score"
      ),
      0
    )
    expect_gt(
      align(
        long,
        short,
        method = method,
        span = "global",
        dissimilarity = "score"
      ),
      0
    )
  }
})

test_that("an affine gap makes score distance larger than identity", {
  a <- "AAAA"
  b <- "AAA"
  d_id <- align(
    a,
    b,
    method = "wfa2",
    match = 0,
    mismatch = 1,
    gap_open = 2,
    gap_extend = 1
  )
  d_sc <- align(
    a,
    b,
    method = "wfa2",
    match = 0,
    mismatch = 1,
    gap_open = 2,
    gap_extend = 1,
    dissimilarity = "score"
  )
  expect_gt(d_sc, d_id)
})

test_that("seq_distmx score mode with cigar uses the score formula", {
  seq <- c(a = "AAAA", b = "AAA")
  mx <- seq_distmx(
    seq,
    threshold = 1,
    dist_config = dist_wfa2(
      match = 0L,
      mismatch = 1L,
      gap_open = 2L,
      gap_extend = 1L,
      dissimilarity = "score"
    ),
    details = "cigar",
    constrain = FALSE,
    parallel_config = parallel_concurrent(1L)
  )
  expect_true("cigar" %in% names(mx))
  expect_equal(nrow(mx), 1L)
  expect_equal(
    mx$dist2,
    align(
      "AAAA",
      "AAA",
      method = "wfa2",
      match = 0,
      mismatch = 1,
      gap_open = 2,
      gap_extend = 1,
      dissimilarity = "score"
    )
  )
  expect_false(identical(
    mx$dist2,
    align(
      "AAAA",
      "AAA",
      method = "wfa2",
      match = 0,
      mismatch = 1,
      gap_open = 2,
      gap_extend = 1
    )
  ))
})

test_that("score mode keeps a cheap-gap pair that identity banding drops", {
  seq <- c(
    a = paste(rep("A", 80L), collapse = ""),
    b = paste(rep("A", 100L), collapse = "")
  )
  cfg_id <- dist_wfa2(
    match = 0L,
    mismatch = 6L,
    gap_open = 0L,
    gap_extend = 1L,
    dissimilarity = "identity"
  )
  cfg_sc <- dist_wfa2(
    match = 0L,
    mismatch = 6L,
    gap_open = 0L,
    gap_extend = 1L,
    dissimilarity = "score"
  )
  hits_id <- seq_distmx(
    seq,
    threshold = 0.1,
    dist_config = cfg_id,
    constrain = TRUE,
    parallel_config = parallel_concurrent(1L)
  )
  hits_sc <- seq_distmx(
    seq,
    threshold = 0.1,
    dist_config = cfg_sc,
    constrain = TRUE,
    parallel_config = parallel_concurrent(1L)
  )
  expect_equal(nrow(hits_id), 0L)
  expect_equal(nrow(hits_sc), 1L)
  expect_lt(hits_sc$dist2[[1]], 0.1)
})
