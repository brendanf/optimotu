library(magrittr)
source("scripts/00_util.R")
source("scripts/01_functions.R")
times <- list()
groups <- c("Penicillium", "Fungi")
for (group in groups) {
  dbseqs <- Biostrings::readDNAStringSet("protaxFungi/addedmodel/sintaxits2train.fa")
  dbseqs <- dbseqs[grepl(group, names(dbseqs))]
  lengths <-  c(1000, 300, 100)
  lengths <- lengths[lengths <= length(dbseqs)]
  # times[[group]] <- tibble::tibble(threads = integer(), minsplit = integer(), length = integer(), time = numeric())
  for (l in lengths) {
    dbseqs <- dbseqs[sample.int(length(dbseqs), l, replace = FALSE)]
    seq_file <- sprintf("temp/%s_ref_%d.fasta", group, l)
    names(dbseqs) <- seq_along(dbseqs) - 1
    Biostrings::writeXStringSet(dbseqs, seq_file)
    hitfile <- sprintf("temp/%s_ref_%d.hits", group, l)
    # usearch_hitlist(dbseqs, threshold = 0.6, ncpu = 1, hits = hitfile,
    #                 usearch = "../GSSP-clust/sh_matching_pub/programs/usearch")
    system2(
      "sh_matching_pub/programs/usearch",
      c(
        "-calc_distmx", seq_file, # input file
        "-tabbedout", hitfile, # output fifo
        "-maxdist", 0.4, # dissimilarity threshold
        "-termdist", 1, # threshold for udist
        "-wordlength", 6,
        "-lopen", "1", # gap opening
        "-lext", "1", # gap extend
        # "-pattern", "111010010111", # pattern gives better result than kmers maybe?
        "-threads", local_cpus()
      )
    )
    times[[group]] <- tidyr::expand_grid(
      threads = c(1, 2, 4, 8)
    ) %>%
      dplyr::mutate(
        length = length(dbseqs),
        version = "RcppThreads",
        time = purrr:::map_dbl(
          threads,
          function(threads) {
            gc()
            system.time(
              single_linkage_matrix_thread(
                file = hitfile,
                seqnames = names(dbseqs),
                dmin = 0.001,
                dmax = 0.4,
                dstep = 0.001,
                threads = threads
              )
            )[3] %T>%
              print()
          }
        )
      ) %>%
       dplyr::bind_rows(
          tibble::tibble(
             threads = 1,
             length = length(dbseqs),
             version = "pool",
             time = {
                gc()
                system.time(
                   single_linkage_pool(
                      file = hitfile,
                      seqnames = names(dbseqs),
                      dmin = 0.001,
                      dmax = 0.4,
                      dstep = 0.001
                   )
                )[3] %T>%
                   print()
             }
          )
       ) %>%
      dplyr::bind_rows(times[[group]], .)
  }
}

library(ggplot2)

times %>%
  dplyr::bind_rows(.id = "group") %>%
  ggplot(aes(x = length, y = time, group = paste(group, version), color = group, shape = version)) +
  geom_point() +
  facet_wrap(~threads, scales = "free_y") +
  scale_x_log10() +
  scale_y_log10()

times %>%
  dplyr::bind_rows(.id = "group") %>%
  # dplyr::mutate(group = factor(group, levels = c("Penicillium", "Trichocomaceae", "Eurotiales", "Eurotiomycetes", "Ascomycota", "Fungi"))) %>%
  ggplot(aes(x = threads, y = time, group = paste(length, version), color = group, shape = version, linetype = version)) +
  geom_point() +
  geom_line() +
  facet_wrap(~group) +
  scale_x_log10() +
  scale_y_log10()

clust1 <- blastclust_repeat(dbseqs, threshold = seq(60, 99.9, 0.1),
                  threshold_name = 600:999, hitlist_method = "usearch",
                  usearch = "protaxFungi/scripts/usearch10.0.240_i86linux32")
clust1 <- lapply(clust1, vapply, trimws, "")
clust1 <- lapply(clust1, lapply, strsplit, " ", fixed = TRUE)
clust1 <- lapply(clust1, lapply, unlist)
clust1 <- lapply(clust1, unname)
clust1 <- lapply(clust1, lapply, as.integer)
clust1 <- lapply(clust1, lapply, sort)
clust1 <- lapply(clust1, function(x) x[order(vapply(x, `[[`, 1L, i = 1))])
clust2 <- lapply(clust2, lapply, as.integer)
clust2 <- lapply(clust2, lapply, sort)
clust2 <- lapply(clust1, function(x) x[order(vapply(x, `[[`, 1L, i = 1))])
clust2 <- rev(clust2)
names(clust2) <- as.character(1000*(1 - as.numeric(names(clust2))))
all.equal(clust1, clust2)

test <- c(rep(TRUE, 20), rep(TRUE, 20))
imin = 1
imax = length(test) + 1;
while (imax > imin + 1) {
  imean = (imin + imax) %/% 2;
  if (test[imean]) {
    imax = imean;
  } else {
    imin = imean;
  }
}

