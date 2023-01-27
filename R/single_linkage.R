#' Title
#'
#' @param file
#' @param seqnames
#' @param method
#' @param thresholds
#' @param dmin
#' @param dmax
#' @param dstep
#' @param threads
#' @param minsplit
#'
#' @return
#' @export
#'
#' @examples
single_linkage = function(
   file,
   seqnames,
   method = c("tree", "matrix"),
   thresholds = NULL,
   dmin = NULL,
   dmax = NULL,
   dstep = NULL,
   threads = 1L,
   minsplit = 1L
) {
   method = match.arg(method)
   if (is.null(thresholds)) {
      if (is.null(dmin) || is.null(dmax) || is.null(dstep)) {
         stop("either 'thresholds' or 'dmin', 'dmax', and 'dstep' must be given")
      }
      if (!is.numeric(dmin) || !is.numeric(dmax) || !is.numeric(dstep)) {
         stop("'dmin', 'dmax', and 'dstep' must all be numbers.")
      }
      if (length(dmin) != 1L || length(dmax) != 1L || length(dstep) != 1L) {
         stop("'dmin', 'dmax', and 'dstep' must all be of length 1.")
      }
      if (is.na(dmin) || is.na(dmax) || is.na(dstep)) {
         stop("'dmin', 'dmax', and 'dstep' may not be NA.")
      }
      if (is.nan(dmin) || is.nan(dmax) || is.nan(dstep)) {
         stop("'dmin', 'dmax', and 'dstep' may not be NaN.")
      }
      if (dstep <= 0) {
         stop("'dstep' must be positive.")
      }
      if (dmax < dmin) {
         stop("'dmax' must be greater than or equal to 'dmin'.")
      }
      switch(
         method,
         tree = single_linkage_pool_uniform(file, seqnames, dmin, dmax, dstep),
         matrix = single_linkage_matrix_uniform(file, seqnames, dmin, dmax, dstep, threads, minsplit)
      )
   } else if (!is.null(dmin) || !is.null(dmax) || !is.null(dstep)) {
      stop("If 'thresholds' is given, 'dmin', 'dmax', and 'dstep' must not be.")
   } else if (!is.numeric(thresholds)) {
      stop("'thresholds' must be numeric.")
   } else if (length(thresholds) == 0) {
      stop("At least one threshold must be given in 'thresholds'.")
   } else if (any(is.na(thresholds))) {
      stop("'thresholds' may not have NA values.")
   } else if (any(is.nan(thresholds))) {
      stop("'thresholds' may not have NaN values.")
   } else if (any(thresholds < 0)) {
      stop("'thresholds' may not contain negative values.")
   } else if (any(utils::head(thresholds, -1) >= utils::tail(thresholds, -1))) {
      stop("'thresholds' must be strictly increasing.")
   } else {
      switch(
         method,
         tree = single_linkage_pool_array(file, seqnames, thresholds),
         matrix = single_linkage_matrix_array(file, seqnames, thresholds, threads, minsplit)
      )
   }
}
