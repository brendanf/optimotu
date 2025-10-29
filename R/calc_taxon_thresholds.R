#' Calculate clustering thresholds for each taxon
#' @param supertaxa (`character` vector) names of taxa at rank `superrank`
#' @param superrank (`character` string) the rank of the taxon within which the
#' sequences are to be clustered
#' @param rank (`character` string) the rank which is to be approximated by
#' clustering
#' @param optima (`data.frame`) optimum clustering thresholds to approximate
#' taxa at different ranks within the taxonomic hierarchy; as returned by
#' `optimize_thresholds()`
#' @return (`numeric`) the optimal threshold for clustering sequences within
#' each `supertaxon`, in order to approximate rank `rank`
#' @keywords internal
superrank_thresholds <- function(supertaxa, superrank, rank, optima) {
  chosen_optima <- optima$rank == rank & optima$superrank == superrank
  chosen_thresholds <- optima$threshold[chosen_optima]
  names(chosen_thresholds) <- optima$supertaxon[chosen_optima]
  unname(chosen_thresholds[supertaxa])
}


#' Calculate clustering thresholds for each taxon, falling back to its ancestor
#' taxa as necessary
#'
#' @param rank (`character` string) the rank within which clustering will be
#' performed
#' @param taxon_table (`data.frame`) "wide" taxonomy table; column "seq_id"
#' gives the sequence ID, and columns {`root_rank()`} to {`tip_rank()`} (e.g.,
#' "kingdom" to "species") give the taxonomy at each rank.
#' @param optima (`data.frame`) optimum clustering thresholds within
#' various taxa;
#' column "rank" gives the rank which is approximated by
#' clustering;
#' "superrank" gives the rank of the taxon within which the clustering threshold
#' was optimized;
#' "supertaxon" gives that taxon name;
#' "threshold" gives the optimum clustering threshold;
#' "value" gives the value of the optimization target at the optimum threshold;
#' optionally, "conf_level" gives a string description of the confidence level;
#' optionally, "metric" gives a string description of the optimization target.
#' Typically, "conf_level" and "measure" are only needed when optimization has
#' been performed for multiple confidence levels or measures, and so the
#' arguments by the same name should also be included.
#' @param ranks (`character` vector) the ranks of the taxon hierarchy, from
#' most inclusive to least inclusive
#' @param conf_level (`character` string or NULL) if given, then only rows of
#' `optima` with this value in the "conf_level" column will be used
#' @param measure (`character` string or NULL) if given, then only rows of
#' `optima` with this value in the "measure" or "metric" column will be used
#' @param default (`character string`) default taxon to define threshold to use
#' when taxonomy is unknown.
#' @param metric (`character` string or NULL; deprecated) if given, then only
#' rows of `optima` with this value in the "measure" or "metric" column will be
#' used
#'
#' @keywords internal
#'
#' @return (named `numeric`, where names are taxa and values are clustering
#' thresholds)
#' unknown sequences, grouped by taxonomy at the parent rank. Sequences where
#' the parent rank is also unknown are in the item named `"_NA_"`.
#' @export
calc_taxon_thresholds <- function(
    rank,
    taxon_table,
    optima,
    ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
    conf_level = NULL,
    measure = NULL,
    default = unique(taxon_table[[ranks[1]]]),
    metric = NULL
) {
  # check arguments
  checkmate::assert_string(rank)
  checkmate::assert_character(ranks, unique = TRUE, any.missing = FALSE)
  checkmate::assert_subset(rank, ranks)
  checkmate::assert_true(rank != ranks[length(ranks)])

  checkmate::assert_data_frame(taxon_table)
  checkmate::assert_names(
    names(taxon_table),
    must.include = c("seq_id", superranks(rank, ranks), rank, subranks(rank, ranks)[1])
  )

  if (nrow(taxon_table) == 0 || all(is.na(taxon_table$seq_id))) {
    return(numeric())
  }
  checkmate::assert_data_frame(optima)
  checkmate::assert_names(
    names(optima),
    must.include = c("rank", "superrank", "supertaxon", "threshold")
  )
  checkmate::assert_string(conf_level, null.ok = TRUE)
  # filter the optima to only include the chosen confidence level
  if (!is.null(conf_level)) {
    checkmate::assert_subset(conf_level, unique(optima$conf_level))
    optima <- optima[optima$conf_level == conf_level, ]
  }

  checkmate::assert_string(measure, null.ok = TRUE)
  checkmate::assert_string(metric, null.ok = TRUE)
  checkmate::assert(checkmate::test_null(metric) || checkmate::test_null(measure),
                    "Only one of 'measure' and 'metric' should be given")
  if (is.null(measure)) {
    measure <- metric
  }
  # filter the optima to only include the chosen measure
  if (!is.null(measure)) {
    if ("metric" %in% names(optima)) {
      checkmate::assert_subset(measure, optima$metric)
      optima <- optima[optima$metric == measure, ]
    } else {
      checkmate::assert_subset(measure, optima$measure)
    optima <- optima[optima$measure == measure, ]
    }
  }
  # ensure the thresholds are represented as fractional (not percent) distances
  # (not similarities)
  optima$threshold <- threshold_as_dist(optima$threshold)

  checkmate::assert_string(default)

  # we only need unique rows, for cases where the taxon is known
  taxon_table <- taxon_table[!is.na(taxon_table[[rank]]), c(superranks(rank, ranks), rank), drop = FALSE]
  taxon_table <- unique(taxon_table)

  # rank is also being used as a superrank (we are clustering inside it)
  my_superranks <- c(superranks(rank, ranks), rank)
  my_child_rank <- subranks(rank, ranks)[1]

  # start by calculating the threshold for the most specific rank
  # and work up to the most general rank

  threshold <- rep(NA_real_, nrow(taxon_table))
  names(threshold) <- taxon_table[[rank]]
  for (r in rev(my_superranks)) {
    sr_thresh <- superrank_thresholds(
      supertaxa = taxon_table[[r]],
      superrank = r,
      rank = my_child_rank,
      optima = optima
    )
    threshold <- ifelse(is.na(threshold), sr_thresh, threshold)
  }

  # if the taxon is unknown, use the default threshold
  na_optima <- optima$rank == my_child_rank &
    optima$supertaxon == default
  na_threshold <- optima$threshold[na_optima]
  names(na_threshold) <- "_NA_"

  c(threshold, na_threshold)
}

#' Convert thresholds which may be distances or similarities to distances
#'
#' @param thresholds (`numeric`) thresholds to convert. These may be distances
#' (0 is identity) with values between 0.0 and 0.5, or between 0.0 and 50.0,
#' in which case they are interpreted as percentages;
#' or similarities (1 or 100 is identity) with values between 0.5 and 1.0, or
#' between 50.0 and 100.0. Note that if *all* percentage distances are < 0.5,
#' they will be misinterpreted as fractional distances.
#' @return (`numeric`) thresholds as distances, in the range from 0.0 to 1.0
#' @export
threshold_as_dist <- function(thresholds) {
  if (any(thresholds > 50, na.rm = TRUE)) {
    1 - 0.01 * thresholds
  } else if (any(thresholds > 1, na.rm = TRUE)) {
    0.01 * thresholds
  } else if (any(thresholds > 0.5, na.rm = TRUE)) {
    1 - thresholds
  } else {
    thresholds
  }
}
