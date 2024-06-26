% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{clust_config}
\alias{clust_config}
\alias{clust_tree}
\alias{clust_matrix}
\alias{clust_index}
\alias{clust_slink}
\title{Configuration for clustering algorithms}
\usage{
clust_config(method = c("tree", "matrix", "index", "slink"), ...)

clust_tree()

clust_matrix(
  binary_search = TRUE,
  fill_method = c("binary", "linear", "topdown")
)

clust_index()

clust_slink()
}
\arguments{
\item{method}{(\code{character} string) The clustering algorithm to use. Options
are \code{"tree"}, \code{"matrix"}, and \code{"index"}}

\item{...}{passed on to variants}

\item{binary_search}{(\code{logical} flag) if \code{TRUE}, use binary search instead
of linear search when determining the currently known minimum distance
between two sequences.  This may be slightly faster when the number of
thresholds is very large.}

\item{fill_method}{(\code{character} string) method to use to determine which
matrix elements must be updated for a sequence.  \code{"linear"} and \code{"binary"}
both update the matrix column with a range memory write (using \code{memcpy})
after using a search (either linear or binary) to determine which range needs
to be updated.  \code{"topdown"} fills element-by element starting at the largest
distance, and is primarily included only as a testing option, since it is
usually slower.}
}
\value{
an object describing the clustering algorithm, to pass to
\code{distmx_cluster()} or \code{seq_cluster()}
}
\description{
Configuration for clustering algorithms
}
\section{Functions}{
\itemize{
\item \code{clust_tree()}: helper function for method \code{"tree"}

\item \code{clust_matrix()}: helper function for method \code{"matrix"}

\item \code{clust_index()}: helper function for method \code{"index"}

\item \code{clust_slink()}: helper function for method \code{"slink"}

}}
\section{Tree method}{
The tree method maintains the state of the clustering algorithm using a
tree structure. When there are more than a few thresholds, it uses less
memory than the matrix methods, and requires fewer operations to update.
However, because its data structure is less cache friendly and it requires
multiple memory accesses to determine whether an incoming pairwise distance
will lead to a cluster update, it is often slower than the matrix-based
methods for small datasets.
}

\section{Matrix method}{
The matrix method maintains the state of the clustering algorithm in a
contiguous matrix structure; for matrix output and when using a single thread
or concurrent parallelism, this does not require any additional memory beyond
what is required to store the output. The matrix method uses the least code,
and despite the higher number of operations required for updates, it is
competitive due to its \eqn{\mathcal{O}(1)}{O(1)} checks on incoming pairwise
distances and good cache locality.
}

\section{Index(ed matrix) method}{
The indexed matrix method uses a linked-list index which maintains a
permutation of the columns of the cluster matrix such that all clusters are
contiguous.  This adds some overhead to updates in order to maintain the
index, but significantly reduces the number of matrix columns which must be
accessed during each update.
}

\section{SLINK method}{
Published by Sibson (1972), the SLINK method is optimally fast and low-memory.
Its only downside is that it can only accept the distances in canonical order;
it is not applicable to methods which may return distances in random order,
such as kmer indexing (e.g., usearch/vsearch) or methods where multiple
asynchronous workers are feeding a single instance of SLINK.
}

\references{
R. Sibson, SLINK: An optimally efficient algorithm for the
single-link cluster method, The Computer Journal, Volume 16, Issue 1, 1973,
Pages 30–34, \doi{10.1093/comjnl/16.1.30}
}
