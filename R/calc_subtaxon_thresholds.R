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
#' @return (named `list` of `double` vectors)
#' @export
calc_subtaxon_thresholds <- function(rank, taxon_table, optima,
                                     ranks = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
                                     conf_level = NULL, measure = NULL,
                                     default = unique(taxon_table[[ranks[1]]]),
                                     metric = NULL) {
  # avoid R CMD check NOTE: no visible binding for global variable
  subrank <- superrank <- supertaxon <- threshold <- NULL

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
    return(list())
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
  my_subranks <- subranks(rank, ranks)

  thresholds <- vector("list", length(my_subranks))
  names(thresholds) <- rev(my_subranks)

  for (child_rank in subranks(rank, ranks)) {
    thresholds[[child_rank]] <- rep(NA_real_, nrow(taxon_table))
    for (r in rev(my_superranks)) {
      sr_thresh <- superrank_thresholds(
        supertaxa = taxon_table[[r]],
        superrank = r,
        rank = child_rank,
        optima = optima
      )
      thresholds[[child_rank]] <- ifelse(
        is.na(thresholds[[child_rank]]),
        sr_thresh,
        thresholds[[child_rank]]
      )
    }
  }

  thresholds <- as.data.frame(thresholds, row.names = taxon_table[[rank]])
  thresholds <- t(thresholds)
  thresholds <- as.data.frame(thresholds)
  thresholds <- unclass(thresholds)
  thresholds <- lapply(thresholds, cummax)
  thresholds <- lapply(thresholds, `names<-`, rev(my_subranks))

  na_optima <- optima$supertaxon == default
  thresholds_na <- optima$threshold[na_optima]
  names(thresholds_na) <- optima$rank[na_optima]
  thresholds_na <- thresholds_na[rev(my_subranks)]
  thresholds_na <- cummax(thresholds_na)
  thresholds[["_NA_"]] <- thresholds_na

  thresholds
}
