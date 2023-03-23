library(magrittr)
source("scripts/00_util.R")
source("scripts/01_functions.R")
times <- list()
groups <- c("Penicillium", "Fungi")
for (group in groups) {
  dbseqs <- Biostrings::readDNAStringSet("protaxFungi/addedmodel/sintaxits2train.fa")
  dbseqs <- dbseqs[grepl(group, names(dbseqs))]
  lengths <-  c(300, 100, 10)
  lengths <- lengths[lengths <= length(dbseqs)]
  # times[[group]] <- tibble::tibble(threads = integer(), minsplit = integer(), length = integer(), time = numeric())
  for (l in lengths) {
    set.seed(1)
    dbseqs <- dbseqs[sample.int(length(dbseqs), l, replace = FALSE)]
    seq_file <- sprintf("temp/%s_ref_%d.fasta", group, l)
    names(dbseqs) <- seq_along(dbseqs) - 1
    if (!file.exists(seq_file)) Biostrings::writeXStringSet(dbseqs, seq_file)
    hitfile <- sprintf("temp/%s_ref_%d.hits", group, l)
    # usearch_hitlist(dbseqs, threshold = 0.6, ncpu = 1, hits = hitfile,
    #                 usearch = "../GSSP-clust/sh_matching_pub/programs/usearch")
    if (!file.exists(hitfile))
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
    hit_conn = file(hitfile, open = "r")
    on.exit(close(hit_conn))
    hitfile2 = paste0(hitfile, "2")
    unlink(hitfile2)
    file.create(hitfile2)
    hit_conn2 = file(hitfile2, open = "a")
    on.exit(close(hit_conn2), add = TRUE)
    dm <- as.dist(matrix(rep(1,l^2), nrow = l))
    cat(group, l)
    while(length(line <- readLines(hit_conn, 1))) {
       distline <- as.numeric(strsplit(line, split = "\t")[[1]])
       if (distline[1] == distline[2]) next
       if (distline[1] < distline[2]) {
          i = distline[1] + 1
          j = distline[2] + 1
       } else {
          j = distline[1] + 1
          i = distline[2] + 1
       }
       dm[l*(i-1) - i*(i-1)/2 + j-i] <- distline[3]
       hc <- hclust(dm, method = "single") %>%
          lapply(0:400 / 1000, cutree, k = NULL, tree = .) %>%
          lapply(function(x) {
             for (i in seq_along(x)) {if (x[i] < i) x[x >= i] = x[x >= i] + 1L}
             x
          }) %>%
          do.call(rbind, .)
       hc <- hc - 1L
      writeLines(line, hit_conn2)
      flush(hit_conn2)
       slmt <- single_linkage_matrix_thread(
          hitfile2,
          names(dbseqs),
          dmin = 0.000,
          dmax = 0.4,
          dstep = 0.001
       )
       slp <- single_linkage_pool(
          hitfile2,
          names(dbseqs),
          dmin = 0.000,
          dmax = 0.4,
          dstep = 0.001
       )
       stopifnot(identical(slmt, hc))
       stopifnot(identical(slmt, slp))
       cat(".")
    }
    cat("\n")
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

