% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/config.R
\name{prealign_config}
\alias{prealign_config}
\alias{prealign_kmer}
\alias{prealign_wfa2}
\alias{prealign_edlib}
\alias{prealign_sneakysnake}
\title{Config helper for prealignment methods}
\usage{
prealign_config(method = c("kmer", "wfa2", "edlib", "sneakysnake"), ...)

prealign_kmer(udist_threshold = 1)

prealign_wfa2()

prealign_edlib()

prealign_sneakysnake()
}
\arguments{
\item{method}{(\code{character} string) prealignment method to use}

\item{...}{passed on to variants}

\item{udist_threshold}{(positive \code{numeric} scalar) kmer distance threshold
for alignment. If <= 1, this is as a fraction of total kmers; if >1 this is
a multiple of the relevant pairwise distance threshold (or 1, whichever is
less). The current implementation uses 8-mers.}
}
\value{
an object representing the prealignment method
}
\description{
Prealignment methods are used to decide which sequence pairs need to have
the full alignment method applied. This is only relevant when the
prealignment method is both faster than the full alignment method, and will
exclude a substantial fraction of all sequence pairs. The most common case
for this is when using linear or affine gap alignment for the actual
distance, but edit distance for the prealignment.  Kmer distance offers
substantial speed increases for the prealignment step, but especially at
dissimilarity thresholds greater than about 0.1, it cannot guarantee that all
relevant sequence pairs are considered.
}
\section{Functions}{
\itemize{
\item \code{prealign_kmer()}: helper function for method \code{"kmer"}

\item \code{prealign_wfa2()}: helper function for method \code{"wfa2"}. As a
prealigner, wfa2 is always configured to use edit distance.

\item \code{prealign_edlib()}: helper function for method \code{"edlib"}

\item \code{prealign_sneakysnake()}: helper function for method \code{"sneakysnake"}

}}
