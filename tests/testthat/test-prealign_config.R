test_that("prealign_config supports kmer, wfa2, and edlib", {
  expect_equal(prealign_config("kmer")$method, "kmer")
  expect_equal(prealign_config("wfa2")$method, "wfa2")
  expect_equal(prealign_config("edlib")$method, "edlib")
})

test_that("prealign_config(\"sneakysnake\") warns and then errors", {
  expect_error(
    expect_warning(
      prealign_config("sneakysnake"),
      "deprecated"
    ),
    "no longer supported"
  )
})
