library(optimotu)
mxfile = "../optimotu/data-raw/cladosporium.mx"
seqnames = as.character(1:6000)
dstep = 0.001
thresholds <- threshold_uniform(0, 0.1, dstep)
microbenchmark::microbenchmark(

  matrix2_bb = distmx_cluster_matrix2_uniform(
    file = mxfile,
    seqnames = seqnames,
    dmin = 0,
    dmax = 0.1,
    dstep = dstep,
    do_binary_search = TRUE,
    fill_method = 2
  ),
  matrix_bb_c2 = optimotu:::distmx_cluster_single(
    file = mxfile,
    seqnames = seqnames,
    threshold_config = thresholds,
    method_config = clust_matrix(TRUE, "binary"),
    parallel_config = parallel_merge(2),
    output_type = "matrix"
  ),
  # matrix_bl_c4 = optimotu:::distmx_cluster_single(
  #   file = mxfile,
  #   seqnames = seqnames,
  #   threshold_config = thresholds,
  #   method_config = clust_matrix(TRUE, "linear"),
  #   parallel_config = parallel_concurrent(4),
  #   output_type = "matrix"
  # ),
  # matrix2_lb = optimotu:::distmx_cluster_single(
  #   file = mxfile,
  #   seqnames = seqnames,
  #   threshold_config = thresholds,
  #   method_config = clust_matrix(FALSE, "binary"),
  #   parallel_config = parallel_concurrent(4),
  #   output_type = "matrix"
  # ),
  # matrix2_ll = optimotu:::distmx_cluster_single(
  #   file = mxfile,
  #   seqnames = seqnames,
  #   threshold_config = thresholds,
  #   method_config = clust_matrix(FALSE, "linear"),
  #   parallel_config = parallel_concurrent(4),
  #   output_type = "matrix"
  # ),
  # matrix2_lL = optimotu:::distmx_cluster_single(
  #   file = mxfile,
  #   seqnames = seqnames,
  #   threshold_config = thresholds,
  #   method_config = clust_matrix(FALSE, "topdown"),
  #   parallel_config = parallel_concurrent(4),
  #   output_type = "matrix"
  # ),
  # pool = distmx_cluster_pool_uniform(
  #   file = mxfile,
  #   seqnames = seqnames,
  #   dmin = 0,
  #   dmax = 0.1,
  #   dstep = dstep,
  #   output_type = "matrix"
  # ),
  # tree = optimotu:::distmx_cluster_single(
  #   file = mxfile,
  #   seqnames = seqnames,
  #   threshold_config = thresholds,
  #   method_config = clust_tree(),
  #   parallel_config = parallel_concurrent(4),
  #   output_type = "matrix"
  # ),
  # index = optimotu:::distmx_cluster_single(
  #   file = mxfile,
  #   seqnames = seqnames,
  #   threshold_config = thresholds,
  #   method_config = clust_index(),
  #   parallel_config = parallel_concurrent(4),
  #   output_type = "matrix"
  # ),
  times = 3,
  check = "equivalent"
)

