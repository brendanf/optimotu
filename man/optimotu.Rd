% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimotu.R, R/optimotu_package.R
\docType{package}
\name{optimotu}
\alias{optimotu-package}
\alias{optimotu}
\title{Taxonomically guided OTU clustering using optimized thresholds}
\usage{
optimotu(
  seqs,
  tax_prob,
  threshold_optima,
  ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
  measure = NULL,
  prob_thresh = NULL,
  dist_config = dist_wfa2(),
  clust_config = clust_slink(),
  parallel_config = parallel_concurrent(1),
  verbose = FALSE
)
}
\arguments{
\item{seqs}{(\code{\link[Biostrings:XStringSet-class]{DNAStringSet}}, file name,
\code{character} vector, or \code{data.frame}) input sequences}

\item{tax_prob}{(\code{data.frame}) taxonomic probabilities for each sequence,
with columns \code{seq_id} giving unique sequence identifiers, \code{rank} giving the
taxonomic rank, \code{taxon} giving the identified taxon name, and, optionally,
\code{prob} giving the probability of the taxon assignment. All sequences in \code{seqs}
must be identified at the first rank listed in \code{ranks} (usually "kingdom");
if this is not the case, it is recommended to add an additional dummy rank,
such as "root_rank", to \code{tax_prob} with the appropriate taxon assignments.
Rows where \code{taxon} or \code{prob} is \code{NA}, or where \code{prob < prob_thresh}, are
considered unknown.}

\item{threshold_optima}{(\code{data.frame}) optimized thresholds, as calculated
by \code{optimize_thresholds()}}

\item{ranks}{(\code{character}) The taxonomic ranks to cluster at, in order from
most inclusive to least inclusive. Defaults to \code{c("kingdom", "phylum", "class", "order", "family", "genus", "species")}. Any extra ranks given in
\code{tax_prob} or \code{threshold_optima} are silently ignored.}

\item{measure}{(\code{character}) The measure to optimize for. Required if
multiple measures are present in \code{threshold_optima}, in which case the given
measure must be included in \code{threshold_optima$measure}.}

\item{prob_thresh}{(\code{numeric}) The probability threshold to use as a cutoff
for accepting taxonomic assignments. Defaults to \code{NULL}, in which case the
all taxonomic assignments in \code{tax_prob} are used. If given, \code{tax_prob} must
include a \code{prob} column.}

\item{dist_config}{(\code{optimotu_dist_config}) configuration for calculating
distances, as returned by \code{dist_config()} or its helpers.}

\item{clust_config}{(\code{optimotu_clust_config}) configuration for clustering,
as returned by \code{clust_config()} or its helpers.}

\item{parallel_config}{(\code{optimotu_parallel_config}) configuration for parallel
processing, as returned by \code{parallel_config()} or its helpers.}

\item{verbose}{(\code{logical} or \code{integer}) print progress messages.}
}
\value{
(\code{data.frame}) with columns \code{seq_id} giving the sequence identifiers
and the taxonomic assignments at each rank in \code{ranks}
}
\description{
Optimize OTU clustering thresholds
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/brendanf/optimotu}
  \item Report bugs at \url{https://github.com/brendanf/optimotu/issues}
}

}
\author{
Brendan Furneaux \href{mailto:brendan.furneaux@gmail.com}{brendan.furneaux@gmail.com}
}
