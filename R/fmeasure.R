#' Calculate multiclass F measure
#'
#' @details 'k' and 'c' can be given in two different representations. In both
#' cases, $n$ objects are clustered, 'k' represents $m$ different "test"
#' clusterings, and 'c' represents the "true" clusters.
#'
#' In the 'matrix' representation, `k` is an integer matrix with $m$ rows and
#' $n$ columns, where `k[i, j]` is the index of the cluster to which element
#' $j$ belongs in the $i$th test clustering. `c` is an integer vector of length
#' $n$, where `c[j]` represents the index of the cluster which element $j$
#' belongs to in the true clustering. Cluster indices need not be dense (i.e.,
#' the presence of clusters with indices 1 and 3 does not imply the presence of
#' a cluster with index 2), and need not match between different test
#' clusterings and the true clustering; all that matters is that two elements
#' with the same value in `k[i,]` (or in `c`) are considered to be clustered
#' together.
#'
#' In the 'set' representation, `k` is a list of lists of integer vectors, where
#' each vector gives the __sorted__ indices of the elements which are included
#' in one cluster.  Then the inner list gives a full clustering of the $n$
#' elements, and the outer list contains $m$ such clusterings. `c` is analogous
#' to one of the inner lists of `k`. The order of clusters does not matter, but
#' the indexing scheme used in `c` and `k` must match; i.e. index 4 must refer to
#' the same element each time it appears in `k` and `c`.  In order to be a
#' coherent clustering, each element should only in exactly one cluster in `c`
#' and each inner list of `k`, but this is not actually checked.
#'
#' The algorithm used for the 'matrix' calculation is considerably faster, and
#' in the future the 'set' representation will probably be converted internally
#' to the 'matrix' representation.
#'
#' @param k (`integer` matrix or `list` of `list`s of `integer` vectors) "test" clusters. See Details
#' @param c (`integer` vector or `list` of `integer` vectors) "true" clusters. See Details
#' @param ncpu (`integer` scalar) number of threads to use
#'
#' @return `numeric` vector giving the F-measure for each row
#' (matrix representation) or top-level element (set representation) of k.
#' @export
fmeasure <- function(k, c, ncpu) {
  if (is.integer(k) && is.matrix(k) && is.integer(c)) {
    fmeasure_matrix(k, c, cpu)
  } else {
    fmeasure_list(k, c, cpu)
  }
}
