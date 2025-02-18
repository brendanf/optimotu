
d <- NULL
test_that("Hamming distance works for short sequences", {
  expect_no_error(
    d <<- seq_search(
      query = c(A="-AAAAACA-GCGT"),
      ref =   c(B="AAAAAACCCGCG-"),
      threshold = 0.4,
      dist_config = dist_hamming(ignore_gaps = TRUE)
    )
  )
  expect_identical(nrow(d), 1L)
  expect_identical(d$seq_id, "A")
  expect_identical(d$ref_id, "B")
  expect_equal(d$dist, 1/10)

  expect_no_error(
    d <<- seq_search(
      query = c(A="-AAAAACA-GCGT"),
      ref =   c(B="AAAAAACCCGCG-"),
      threshold = 0.4,
      dist_config = dist_hamming(ignore_gaps = FALSE)
    )
  )
  expect_identical(nrow(d), 1L)
  expect_identical(d$seq_id, "A")
  expect_identical(d$ref_id, "B")
  expect_equal(d$dist, 2/11)

})


test_that("Hamming distance works for long sequences", {
  expect_no_error(
    d <<- seq_search(
      query = c(A="--ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAAACA-GCGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
      ref   = c(B="GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAAACCCGCGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT--"),
      threshold = 0.4,
      dist_config = dist_hamming(ignore_gaps = TRUE)
    )
  )
  expect_identical(nrow(d), 1L)
  expect_identical(d$seq_id, "A")
  expect_identical(d$ref_id, "B")
  expect_equal(d$dist, 1/138)

  expect_no_error(
    d <<- seq_search(
      query = c(A="--ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAAACA-GCGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"),
      ref =   c(B="CTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAAACCCGCGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT--"),
      threshold = 0.4,
      dist_config = dist_hamming(ignore_gaps = FALSE)
    )
  )
  expect_identical(nrow(d), 1L)
  expect_identical(d$seq_id, "A")
  expect_identical(d$ref_id, "B")
  expect_equal(d$dist, 2/139)

})
