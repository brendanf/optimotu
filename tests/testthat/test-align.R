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
