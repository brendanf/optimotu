test_that("align gives proper values for global and extension", {
  expect_equal(align("AAACGCGTCG", "AAACGCGT", method = "wfa2", span = "global"), 0.2)
  expect_equal(align("AAACGCGTCG", "AAACGCGT", method = "wfa2", span = "extend"), 0)
  expect_equal(align("AAACGCGTCG", "AAACGCGT", method = "edlib", span = "global"), 0.2)
  expect_equal(align("AAACGCGT", "AAACGCGTCG", method = "edlib", span = "extend"), 0)
})

