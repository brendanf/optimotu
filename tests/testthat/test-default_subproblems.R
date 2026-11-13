testthat::test_that("default_subproblems fast methods use square nearest 4x", {
  testthat::expect_equal(default_subproblems(1L, "hamming"), 1L)
  testthat::expect_equal(default_subproblems(16L, "hamming"), 64L)
  testthat::expect_equal(default_subproblems(15L, "hamming"), 64L)
  testthat::expect_equal(default_subproblems(30L, "hamming"), 121L)
  testthat::expect_equal(default_subproblems(12L, "usearch"), 49L)
  testthat::expect_equal(default_subproblems(16L, "aa"), 64L)
  testthat::expect_equal(default_subproblems(16L, "evo"), 64L)
})

testthat::test_that("default_subproblems slow methods maximize last-wave fill", {
  testthat::expect_equal(default_subproblems(16L, "ksw2"), 16L)
  testthat::expect_equal(default_subproblems(15L, "ksw2"), 25L)
  testthat::expect_equal(default_subproblems(30L, "wfa2"), 81L)
  testthat::expect_equal(default_subproblems(20L, "edlib"), 36L)
  testthat::expect_equal(default_subproblems(12L, "hybrid"), 36L)
})

testthat::test_that("resolve_subproblems keeps explicit values and floors at threads", {
  testthat::expect_equal(
    optimotu:::resolve_subproblems(16L, "hamming", NULL),
    64L
  )
  testthat::expect_equal(
    optimotu:::resolve_subproblems(16L, "hamming", 36L),
    36L
  )
  testthat::expect_equal(
    optimotu:::resolve_subproblems(16L, "hamming", 4L),
    16L
  )
})

testthat::test_that("parallel_* store subproblems only when set", {
  pc <- parallel_merge(4L)
  testthat::expect_null(pc$subproblems)
  pc2 <- parallel_merge(4L, subproblems = 36L)
  testthat::expect_equal(pc2$subproblems, 36L)
  pc3 <- parallel_concurrent(8L, subproblems = 25L)
  testthat::expect_equal(pc3$subproblems, 25L)
})
