test_that("submodule pins match the gitlinks", {
  root <- normalizePath(testthat::test_path("../.."), mustWork = TRUE)
  if (!file.exists(file.path(root, ".git"))) {
    skip("not a git checkout")
  }
  script <- file.path(root, "tools", "fetch-submodules.sh")
  out <- system2("sh", c(script, "--check"), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (status != 0L) {
    message(paste(out, collapse = "\n"))
  }
  expect_identical(status, 0L)
})
