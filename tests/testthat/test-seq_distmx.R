test_that("seq_distmx_internal works", {
  testseqs <- seqs <-
    c(
      "ACGT",
      "ACGA",
      "ACGG",
      "ACGC",
      "ACGT",
      "AGCA"
    )
  names(testseqs) <- paste0("seq", seq_along(seqs))
  for (method in c("usearch", "wfa2", "hamming", "edlib")) {
    for (threshold in c(0.1, 0.3, 0.55)) {
      for (threads in c(1, 4)) {
        for (detail in c("cigar", "gapstats", "none")) {
          if (method == "hamming" && detail == "cigar") next
          # cat("\nmethod =", method, ", threshold =", threshold, ", threads =", threads, ", detail =", detail, "\n")

          distmx <- seq_distmx(
            seq = testseqs,
            dist_config = optimotu::dist_config(method = method),
            threshold = threshold,
            parallel_config = optimotu::parallel_concurrent(threads),
            detail = detail,
            # these parameters are for usearch, to make it work in this test
            # normally it cannot be run with such short sequences
            fulldp = "",
            gapopen = "'*E'"
          )
          expect_true(is.data.frame(distmx))
          if (threshold == 0.1) {
            expect_equal(nrow(distmx), 1)
          } else if (threshold == 0.3) {
            expect_equal(nrow(distmx), 10)
          } else if (threshold == 0.5) {
            expect_equal(nrow(distmx), 15)
          }
          if (detail == 2) {
            expect_type(distmx$cigar, "character")
          } else if (detail == 1) {
            expect_type(distmx$align_length, "integer")
            expect_type(distmx$n_insert, "integer")
            expect_type(distmx$n_delete, "integer")
            expect_type(distmx$max_insert, "integer")
            expect_type(distmx$max_delete, "integer")
          }
        }
      }
    }
  }
})
