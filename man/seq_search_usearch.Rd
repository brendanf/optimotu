% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/seq_search_usearch.R
\name{seq_search_usearch}
\alias{seq_search_usearch}
\title{Search for the best match(es) to query sequences in a set of reference
sequences using USEARCH}
\usage{
seq_search_usearch(
  query,
  ref,
  threshold,
  query_id = NULL,
  ref_id = NULL,
  parallel_config = parallel_concurrent(1),
  usearch = find_usearch(),
  verbose = FALSE,
  ...
)
}
\arguments{
\item{query}{(\code{character}, \code{data.frame}, or \code{XStringSet}) sequences to search}

\item{ref}{(\code{character}, \code{data.frame}, or \code{XStringSet}) reference sequences to
search against}

\item{threshold}{(\code{numeric} scalar) maximum distance to consider a match, in
[0, 1] where 0 is identical.}

\item{query_id}{(\code{character} vector) names for the query sequences.  If they
are already named, this will replace the names.  Has no effect if \code{query} is
a filename.}

\item{ref_id}{(\code{character} vector) names for the reference sequences.  If
they are already named, this will replace the names.  Has no effect if \code{seq}
is a filename.}

\item{parallel_config}{(\code{optimotu_parallel_config}) configuration for parallel
processing, as returned by \code{parallel_config()} or its helpers.}

\item{usearch}{(\code{character} scalar) path to USEARCH executable}

\item{verbose}{(\code{logical} or \code{integer} scalar) print progress messages.}

\item{...}{passed to methods (currently inactive)}
}
\value{
(\code{data.frame}) with columns "seq_id" (\code{character}),
"ref_id" (\code{character}), and "dist" (\code{numeric}) giving the distance between
the query and reference.
}
\description{
This function uses a unix pipe to read the output of the
USEARCH "\code{usearch_global}" command. USEARCH (version 8.0' or higher) should
be installed separately; it is available with a free license for most users
at https://www.drive5.com/usearch/.
}
