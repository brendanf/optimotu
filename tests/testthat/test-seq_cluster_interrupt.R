simulate_interrupt_sequences <- function(
    n_clusters = 8L,
    members_per_cluster = 15L,
    seq_len = 260L,
    centroid_sub_frac = 0.10,
    member_sub_frac = 0.03
) {
  alphabet <- c("A", "C", "G", "T")

  mutate_substitutions <- function(seq, frac) {
    out <- strsplit(seq, split = "", fixed = TRUE)[[1]]
    n_mut <- max(1L, as.integer(round(length(out) * frac)))
    idx <- sample.int(length(out), n_mut, replace = FALSE)
    for (i in idx) {
      out[i] <- sample(alphabet[alphabet != out[i]], size = 1L)
    }
    paste0(out, collapse = "")
  }

  ancestor <- paste0(sample(alphabet, seq_len, replace = TRUE), collapse = "")
  centroids <- vapply(
    seq_len(n_clusters),
    function(i) mutate_substitutions(ancestor, centroid_sub_frac),
    character(1)
  )
  seq <- unlist(lapply(seq_len(n_clusters), function(cluster_id) {
    vapply(
      seq_len(members_per_cluster),
      function(member_id) {
        mutate_substitutions(centroids[[cluster_id]], member_sub_frac)
      },
      character(1)
    )
  }), use.names = FALSE)
  names(seq) <- as.character(seq_len(length(seq)))
  seq
}

testthat::test_that("seq_cluster handles SIGINT without crashing process", {
  testthat::skip_on_cran()
  testthat::skip_if(.Platform$OS.type != "unix", "requires unix signals")
  testthat::skip_if_not(
    identical(tolower(Sys.getenv("OPTIMOTU_RUN_INTERRUPT_TESTS")), "true"),
    "set OPTIMOTU_RUN_INTERRUPT_TESTS=true to run interrupt test"
  )

  set.seed(12)
  seq <- simulate_interrupt_sequences()
  job <- parallel::mcparallel({
    tryCatch(
      {
        seq_cluster(
          seq = seq,
          dist_config = dist_hamming(),
          threshold_config = threshold_uniform(0, 0.2, 0.005),
          clust_config = clust_tree(),
          parallel_config = parallel_concurrent(1L),
          output_type = "matrix",
          verbose = FALSE
        )
        "completed"
      },
      interrupt = function(e) "interrupted",
      error = function(e) conditionMessage(e)
    )
  }, silent = TRUE)

  Sys.sleep(0.2)
  tools::pskill(job$pid, tools::SIGINT)
  res <- parallel::mccollect(job, wait = TRUE, timeout = 60)

  testthat::expect_true(!is.null(res))
  testthat::expect_true(
    identical(res[[1]], "interrupted") || identical(res[[1]], "completed"),
    info = sprintf("child result=%s", as.character(res[[1]]))
  )
})
