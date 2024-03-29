% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{parallel_config}
\alias{parallel_config}
\alias{parallel_merge}
\alias{parallel_concurrent}
\alias{parallel_hierarchical}
\title{Configuration for parallelization options}
\usage{
parallel_config(method = c("merge", "concurrent", "hierarchical"), ...)

parallel_merge(threads)

parallel_concurrent(threads)

parallel_hierarchical(threads, shards)
}
\arguments{
\item{method}{(\code{character} string) parallelization method.}

\item{...}{passed on to variants}

\item{threads}{(positive \code{integer} scalar) total number of threads to use}

\item{shards}{(positive \code{integer} scalar) number of independent working
units for clustering. Must be less than or equal to \code{threads} (and in order
for the result to actually be hierarchical, should be at least 2 and at most
\code{threads}/2)}
}
\value{
an object representing the thresholds
}
\description{
Configuration for parallelization options
}
\section{Functions}{
\itemize{
\item \code{parallel_merge()}: helper function for method \code{"merge"}

\item \code{parallel_concurrent()}: helper function for method \code{"concurrent"}

\item \code{parallel_hierarchical()}: helper function for method \code{"hierarchical"}

}}
\section{Merge method}{
In the merge method, each thread works on an independent clustering of the
data, based on its own (disjoint) subset of the distance matrix. When each
thread finishes clustering, it merges its results into the master clustering.
This avoids concurrency collisions between the threads during the main
clustering, although collisions can still occur if multiple threads try to
merge at the same time.
}

\section{Concurrent method}{
In the concurrent method, all threads work jointly on the same clustering of
the data. Although many threads can simultaneously read the clustering to
determine whether a new pairwise distance will lead to an update (most do
not), only one thread can update the clustering at a time, so this method
can lead to more concurrency collisions when many threads are in use.
However, sharing the state of the clustering between threads leads to fewer
total updates than the merge method.
}

\section{Hierarchical method}{
The hierarchical method is a combination of several "shards", each of which
in turn has multiple threads. Each shard has one clustering object, which the
threads within the shard update concurrently. When all threads in the shard
have finished, then the results from the shard are merged into the master
clustering. This is probably the most efficient method when there are very
many threads, but the optimal number of shards varies between data sets.
}

