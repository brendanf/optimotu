test_that("clustering works on a tiny dataset", {
  thresholds <- c(species = 0.017, genus = 0.067, family = 0.084, order = 0.148)
  for (i in 1:10)
    expect_equal(
      distmx_cluster(
        distmx = test_path("fifo3993281b53243b"),
        names = c("ASV59674", "ASV60755"),
        threshold_config = threshold_set(thresholds),
        clust_config = clust_tree(),
        parallel_config = parallel_concurrent(2)
      ),
      matrix(
        c(0, 0, 0, 0, 1, 1, 1, 0),
        ncol = 2,
        dimnames = list(
          c("species", "genus", "family", "order"),
          NULL
        )
      )
    )
})
