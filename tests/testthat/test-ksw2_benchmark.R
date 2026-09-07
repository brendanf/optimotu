# Opt-in KSW2 vs WFA2 performance benchmark on real ITS ASV data.
#
# Thresholds are intentionally large (0.4) to force many alignments and stress
# banded DP. Calibration grows n from 2000 until single-thread seq_distmx takes
# at least 5 minutes (or the FASTA is exhausted).
#
# Run with:
#   OPTIMOTU_RUN_KSW2_BENCHMARK=true Rscript -e 'devtools::test(filter="ksw2_benchmark")'

testthat::test_that("ksw2 vs wfa2 seq_distmx/seq_cluster thread scaling", {
  testthat::skip_on_cran()
  testthat::skip_if_not(
    identical(tolower(Sys.getenv("OPTIMOTU_RUN_KSW2_BENCHMARK")), "true"),
    "set OPTIMOTU_RUN_KSW2_BENCHMARK=true to run ksw2 benchmark"
  )
  testthat::skip_if_not_installed("microbenchmark")
  testthat::skip_if_not_installed("Biostrings")

  fasta <- file.path("..", "lifeplan_ITS", "output", "all_asv.fasta.gz")
  # Resolve relative to package root whether tests run from pkg or tests/
  if (!file.exists(fasta)) {
    fasta <- file.path(
      "..",
      "..",
      "lifeplan_ITS",
      "output",
      "all_asv.fasta.gz"
    )
  }
  if (!file.exists(fasta)) {
    # Prefer absolute path next to this package checkout
    pkg_root <- testthat::test_path("..", "..")
    fasta <- file.path(
      dirname(normalizePath(pkg_root)),
      "lifeplan_ITS",
      "output",
      "all_asv.fasta.gz"
    )
  }
  testthat::skip_if_not(
    file.exists(fasta),
    paste0("missing benchmark FASTA: ", fasta)
  )

  threshold <- 0.4
  target_secs <- 5 * 60
  max_threads <- min(20L, parallel::detectCores(logical = TRUE))
  if (is.na(max_threads) || max_threads < 1L) {
    max_threads <- 1L
  }
  threads <- c(1L, 2L, 4L, 8L, 12L, 16L, 20L)
  threads <- threads[threads <= max_threads]


  # Calibration: grow n until 1-thread WFA2 seq_distmx >= 5 minutes
  n <- 500L
  calibrated <- FALSE
  seqs <- NULL
  while (!calibrated) {
    message(sprintf("Calibrating with n=%d ...", n))
    dna <- Biostrings::readDNAStringSet(fasta, nrec = n)
    testthat::skip_if(
      length(dna) < 2L,
      "benchmark FASTA has fewer than 2 sequences"
    )
    if (length(dna) < n) {
      message(sprintf(
        "FASTA exhausted at %d sequences; using all available",
        length(dna)
      ))
      n <- length(dna)
      calibrated <- TRUE
    }
    seqs <- as.character(dna)
    names(seqs) <- names(dna)
    if (is.null(names(seqs)) || anyDuplicated(names(seqs))) {
      names(seqs) <- as.character(seq_along(seqs))
    }

    t0 <- proc.time()[["elapsed"]]
    invisible(seq_distmx(
      seq = seqs,
      threshold = threshold,
      dist_config = dist_wfa2(mismatch = 2, gap_open = 3, gap_extend = 1),
      parallel_config = parallel_concurrent(1L),
      details = "none",
      constrain = TRUE
    ))
    elapsed <- proc.time()[["elapsed"]] - t0
    message(sprintf("n=%d single-thread seq_distmx (wfa2): %.1f s", n, elapsed))
    if (elapsed >= target_secs || calibrated) {
      calibrated <- TRUE
    } else {
      n <- as.integer(ceiling(n * 1.5))
    }
  }
  message(sprintf("Using n=%d sequences for benchmark", length(seqs)))

  thresh_cfg <- threshold_uniform(0, 0.4, 0.001)
  backends <- list(
    wfa2 = dist_wfa2(mismatch = 2, gap_open = 3, gap_extend = 1),
    ksw2 = dist_ksw2(mismatch = 2, gap_open = 3, gap_extend = 1)
  )
  rows <- list()

  time_once <- function(expr) {
    # One timed replicate; microbenchmark adds overhead for long runs
    t0 <- proc.time()[["elapsed"]]
    force(expr)
    proc.time()[["elapsed"]] - t0
  }

  for (backend in names(backends)) {
    for (thr in threads) {
      message(sprintf("seq_distmx %s threads=%d", backend, thr))
      secs <- time_once({
        invisible(seq_distmx(
          seq = seqs,
          threshold = threshold,
          dist_config = backends[[backend]],
          parallel_config = parallel_concurrent(thr),
          details = "none",
          constrain = TRUE
        ))
      })
      rows[[length(rows) + 1L]] <- data.frame(
        backend = backend,
        task = "seq_distmx",
        threads = thr,
        seconds = secs,
        n = length(seqs),
        stringsAsFactors = FALSE
      )

      message(sprintf("seq_cluster %s threads=%d", backend, thr))
      secs <- time_once({
        invisible(seq_cluster(
          seq = seqs,
          dist_config = backends[[backend]],
          threshold_config = thresh_cfg,
          clust_config = clust_tree(),
          parallel_config = parallel_concurrent(thr),
          output_type = "matrix",
          which = TRUE,
          verbose = FALSE
        ))
      })
      rows[[length(rows) + 1L]] <- data.frame(
        backend = backend,
        task = "seq_cluster",
        threads = thr,
        seconds = secs,
        n = length(seqs),
        stringsAsFactors = FALSE
      )
    }
  }

  result <- do.call(rbind, rows)
  print(result)

  out_dir <- testthat::test_path("_benchmarks")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, "ksw2_vs_wfa2.csv")
  utils::write.csv(result, out_file, row.names = FALSE)
  message("Wrote ", out_file)

  testthat::expect_true(nrow(result) > 0)
  testthat::expect_true(all(result$seconds > 0))
})
