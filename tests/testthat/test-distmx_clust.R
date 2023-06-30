n <- NULL
subsets <- NULL
distmx <- NULL
subsetdistmx <- NULL
dist_table <- NULL
testthat::test_that("dist matrix generation works", {
  testthat::expect_no_error({
    set.seed(1)
    k <- 5L
    m <- 5L
    n <<- k*m
    # generate centroids for 5 clusters in 2d space
    centroids = data.frame(
      x = runif(k),
      y = runif(k)
    )
    # for each cluster, generate 5 random points around it
    points <- lapply(centroids, rep, each = m) |>
      lapply(rnorm, n = n, sd = 0.1) |>
      as.data.frame()

    # take some subsets for testing subset clustering
    subsets <<- lapply(5:20, sample.int, n = n)

    # euclidian distance matrix around the points, rounded to precision of 0.05
    distmx <<- round(dist(points) / 0.05) * 0.05
    subsetdistmx <<- lapply(subsets, \(i) round(dist(as.matrix(points)[i,]) / 0.05)*0.05)


    # convert distance matrix to "sparse" format
    dist_table <<- data.frame(
      seq1 = rep(0:(n-2), (n-1):1),
      seq2 = unlist(lapply(1:(n-1), seq, to = (n-1))),
      dist = unclass(distmx)
    )
    dist_table <<- dist_table[dist_table$dist <= 1,]
    # add tautological 0 distances (USEARCH produces these)
    dist_table <<- rbind(
      data.frame(
        seq1 = 1L:n - 1L,
        seq2 = 1L:n - 1L,
        dist = 0
      ),
      dist_table
    )

    # create 5 shuffles of the dist matrix
    dist_table <<- c(
      list(dist_table),
      replicate(5, dist_table[sample(nrow(dist_table)),], simplify = FALSE)
    )
  })
})

# write as temp files
distmx_file <- character(length(dist_table))
for (i in seq_along(dist_table)) {
  distmx_file[i] <- tempfile(fileext = ".distmx")
  testthat::test_that(sprintf("writing file %s works", distmx_file[i]), {
    testthat::expect_no_error({
      write.table(
        x = dist_table[i],
        file = distmx_file[i],
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE
      )
    })
  })
}
# withr::defer(unlink(distmx_file), teardown_env())

hclust2matrix <- function(distmx, thresholds) {
  out <- hclust(pmax(distmx - 1e-8, 0), method = "single") |>
    lapply(thresholds, cutree, k = NULL, tree = _) |>
    lapply(function(x) {
      for (i in seq_along(x)) {if (x[i] < i) x[x >= i] = x[x >= i] + 1L}
      x
    }) |>
    do.call(rbind, args = _)
  out - 1L
}

hclust_matrix <- NULL
testthat::test_that("hclust works for test matrix", {
  testthat::expect_no_error(
    hclust_matrix <<- hclust2matrix(distmx, 0:20 / 20)
  )
})
subset_hclust_matrix <- NULL
testthat::test_that("hclust works for test subset matrices", {
  testthat::expect_no_error(
    subset_hclust_matrix <<- lapply(subsetdistmx, hclust2matrix, 0:20 / 20)
  )
})

thresholds <- NULL
testthat::test_that("threshold definitions work", {
  testthat::expect_no_error({
    tn <- as.character(0:20/20)

    thresholds <<- list(
      set_thresh = threshold_set(0:20 / 20),
      uniform_thresh = threshold_uniform(0, 1, 0.05),
      lookup_thresh = threshold_lookup(0:20 / 20, 0.05),
      set_thresh_named = threshold_set(0:20 / 20, tn),
      uniform_thresh_named = threshold_uniform(0, 1, 0.05, tn),
      lookup_thresh_named = threshold_lookup(0:20 / 20, 0.05, tn)
    )
  })
})

algorithms <- NULL
testthat::test_that("algorithm definitions work", {
  testthat::expect_no_error({
    algorithms <<- list(
      tree = clust_tree(),
      matrix_binary_binary = clust_matrix(binary_search = TRUE, fill_method = "binary"),
      matrix_binary_linear = clust_matrix(binary_search = TRUE, fill_method = "linear"),
      matrix_binary_topdown = clust_matrix(binary_search = TRUE, fill_method = "topdown"),
      matrix_linear_binary = clust_matrix(binary_search = FALSE, fill_method = "binary"),
      matrix_linear_linear = clust_matrix(binary_search = FALSE, fill_method = "linear"),
      matrix_linear_topdown = clust_matrix(binary_search = FALSE, fill_method = "topdown"),
      index = clust_index()
    )
  })
})

parallels <- NULL
testthat::test_that("parallelism definitions work", {
  testthat::expect_no_error({
    parallels <<- list(
      serial = parallel_concurrent(1),
      concurrent = parallel_concurrent(4),
      merge = parallel_merge(4),
      hierarchical = parallel_hierarchical(2, 2)
    )
  })
})

for (p in names(parallels)) {
  for (t in names(thresholds)) {
    for (a in names(algorithms)) {
      for (i in seq_along(distmx_file)) {
        # cat(sprintf(
        #   "%s distmx_cluster_single %s method with %s thresholds agrees with hclust\n",
        #   p, a, t
        # ))
        test_that(
          sprintf(
            "%s distmx_cluster_single %s method with %s thresholds agrees with hclust (permutation %i)",
            p, a, t, i
          ),
          {
            expect_equal(
              distmx_cluster(
                distmx = distmx_file[i],
                names = as.character(1:n),
                threshold_config = thresholds[[t]],
                clust_config = algorithms[[a]],
                parallel_config = parallels[[p]],
                output_type = "matrix"
              ),
              hclust_matrix,
              ignore_attr = TRUE
            )
          }
        )
        # cat(
        #   sprintf(
        #     "%s distmx_cluster_multi %s method with %s thresholds agrees with hclust\n",
        #     p, a, t
        #   )
        # )
        test_that(
          sprintf(
            "%s distmx_cluster_multi %s method with %s thresholds agrees with hclust (permutation %i)",
            p, a, t, i
          ),
          {
            expect_equal(
              distmx_cluster(
                distmx = distmx_file[i],
                names = as.character(1:n),
                which = lapply(subsets, as.character),
                threshold_config = thresholds[[t]],
                clust_config = algorithms[[a]],
                parallel_config = parallels[[p]],
                output_type = "matrix"
              ),
              subset_hclust_matrix,
              ignore_attr = TRUE
            )
          }
        )
      }
    }
  }
}
