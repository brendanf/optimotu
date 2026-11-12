d <- NULL
test_that("Hamming distance works for short sequences", {
  expect_no_error(
    d <<- seq_search(
      query = c(A = "-AAAAACA-GCGT"),
      ref = c(B = "AAAAAACCCGCG-"),
      threshold = 0.4,
      dist_config = dist_hamming(ignore_gaps = TRUE)
    )
  )
  expect_identical(nrow(d), 1L)
  expect_identical(d$seq_id, "A")
  expect_identical(d$ref_id, "B")
  expect_equal(d$dist, 1 / 10)

  expect_no_error(
    d <<- seq_search(
      query = c(A = "-AAAAACA-GCGT"),
      ref = c(B = "AAAAAACCCGCG-"),
      threshold = 0.4,
      dist_config = dist_hamming(ignore_gaps = FALSE)
    )
  )
  expect_identical(nrow(d), 1L)
  expect_identical(d$seq_id, "A")
  expect_identical(d$ref_id, "B")
  expect_equal(d$dist, 2 / 11)
})


test_that("Hamming distance works for long sequences", {
  expect_no_error(
    d <<- seq_search(
      query = c(
        A = "--ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAAACA-GCGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
      ),
      ref = c(
        B = "GTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAAACCCGCGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT--"
      ),
      threshold = 0.4,
      dist_config = dist_hamming(ignore_gaps = TRUE)
    )
  )
  expect_identical(nrow(d), 1L)
  expect_identical(d$seq_id, "A")
  expect_identical(d$ref_id, "B")
  expect_equal(d$dist, 1 / 138)

  expect_no_error(
    d <<- seq_search(
      query = c(
        A = "--ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAAACA-GCGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC"
      ),
      ref = c(
        B = "CTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAAAAACCCGCGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT--"
      ),
      threshold = 0.4,
      dist_config = dist_hamming(ignore_gaps = FALSE)
    )
  )
  expect_identical(nrow(d), 1L)
  expect_identical(d$seq_id, "A")
  expect_identical(d$ref_id, "B")
  expect_equal(d$dist, 2 / 139)
})

test_that("ignore_gaps FALSE keeps a full word at 64-base ends", {
  for (n in c(64L, 128L, 192L)) {
    same <- strrep("A", n)
    hit <- seq_search(
      query = c(q = same),
      ref = c(r = same),
      threshold = 0.5,
      dist_config = dist_hamming(ignore_gaps = FALSE)
    )
    expect_equal(hit$dist, 0, info = n)

    mismatch <- paste0(strrep("A", n - 1L), "C")
    missed <- seq_search(
      query = c(q = same),
      ref = c(r = mismatch),
      threshold = 0.5,
      dist_config = dist_hamming(ignore_gaps = FALSE)
    )
    expect_equal(missed$dist, 1 / n, info = n)
  }

  # Partial end mask, and ignore_gaps TRUE at a 64-base end.
  partial <- strrep("A", 65L)
  partial_hit <- seq_search(
    query = c(q = partial),
    ref = c(r = partial),
    threshold = 0.5,
    dist_config = dist_hamming(ignore_gaps = FALSE)
  )
  expect_equal(partial_hit$dist, 0)
  full_word <- strrep("A", 64L)
  gaps_ignored <- seq_search(
    query = c(q = full_word),
    ref = c(r = full_word),
    threshold = 0.5,
    dist_config = dist_hamming(ignore_gaps = TRUE)
  )
  expect_equal(gaps_ignored$dist, 0)
})

test_that("Hamming search packing skips lowercase without overflowing", {
  # Packed width is non-lowercase columns. Lowercase is skipped as an
  # insertion, matching clustering. Equal packed width must not over-read.
  q <- c(q1 = "ACGTACGT")
  r <- c(r1 = "ACgtGTACGT")
  d <- seq_search(
    query = q,
    ref = r,
    threshold = 0.5,
    dist_config = dist_hamming(),
    parallel_config = parallel_concurrent(1)
  )
  expect_identical(nrow(d), 1L)
  expect_equal(d$dist, 0)

  expect_error(
    seq_search(
      query = c(q1 = "ACGTACGT"),
      ref = c(r1 = "ACgtGT"),
      threshold = 0.5,
      dist_config = dist_hamming()
    ),
    "same length"
  )
})
