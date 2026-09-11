test_that("align gives proper values for global and extension", {
  expect_equal(
    align("AAACGCGTCG", "AAACGCGT", method = "wfa2", span = "global"),
    0.2
  )
  expect_equal(
    align("AAACGCGTCG", "AAACGCGT", method = "wfa2", span = "extend"),
    0
  )
  expect_equal(
    align("AAACGCGTCG", "AAACGCGT", method = "edlib", span = "global"),
    0.2
  )
  expect_equal(
    align("AAACGCGT", "AAACGCGTCG", method = "edlib", span = "extend"),
    0
  )
  expect_equal(
    align("AAACGCGTCG", "AAACGCGT", method = "ksw2", span = "global"),
    0.2
  )
  expect_equal(
    align("AAACGCGTCG", "AAACGCGT", method = "ksw2", span = "extend"),
    0
  )
})

test_that("wfa2 matches edlib when the first mismatch is past 4 bytes", {
  # WFA2 packed extend XORs 8-byte blocks and uses trailing-zero count.
  # On Windows, __builtin_ctzl() is 32-bit, so a first mismatch in
  # bytes 5-8 of a block stopped match-extend after 4 bases.
  a <- paste(c(rep("A", 5L), "C", rep("G", 10L)), collapse = "")
  b <- paste(c(rep("A", 5L), "T", rep("G", 10L)), collapse = "")
  expect_equal(nchar(a), 16L)
  expect_equal(
    align(a, b, method = "wfa2", span = "global"),
    1 / 16
  )
  expect_equal(
    align(a, b, method = "wfa2", span = "global"),
    align(a, b, method = "edlib", span = "global")
  )
})

test_that("ksw2 agrees with edlib under edit-distance defaults", {
  a <- "AAACGCGTCG"
  b <- "AAACGCGT"
  expect_equal(
    align(a, b, method = "ksw2", span = "global"),
    align(a, b, method = "edlib", span = "global")
  )
  expect_equal(
    align("ACGT", "ACGA", method = "ksw2"),
    align("ACGT", "ACGA", method = "edlib")
  )
})
