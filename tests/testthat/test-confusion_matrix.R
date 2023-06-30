x <- NULL
testthat::test_that("creating matrix x succeeds", {
  expect_no_error(
    x <<- matrix(
      c(1L, 2L, 2L, 4L, 4L,
        1L, 2L, 2L, 1L, 2L),
      nrow = 2,
      byrow = TRUE
    )
  )
})

y <- NULL
testthat::test_that("creating matrix y succeeds", {
  expect_no_error(
    y <<- c(1, 1, 1, 4, 1)
  )
})

cm <- NULL
testthat::test_that("confusion_matrix() works", {
  testthat::expect_equal(
    cm <<-confusion_matrix(x, y, 1),
    data.frame(
      TP = c(1L, 3L),
      FP = c(1L, 1L),
      FN = c(5L, 3L),
      TN = c(3L, 3L)
    )
  )
})

testthat::test_that("confusion matrix based metrics work", {
  testthat::expect_equal(matthews_correlation_coefficient(x, y, 1), matthews_correlation_coefficient(cm))
  testthat::expect_equal(rand_index(x, y, 1), rand_index(cm))
  testthat::expect_equal(adjusted_rand_index(x, y, 1), adjusted_rand_index(cm))
  testthat::expect_equal(fowlkes_mallow_index(x, y, 1), fowlkes_mallow_index(cm))
})
