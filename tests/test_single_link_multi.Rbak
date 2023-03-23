library(magrittr)
source("scripts/00_util.R")
source("scripts/01_functions.R")
times <- list()
groups <- c("Penicillium", "Fungi")
for (group in groups) {
  dbseqs <- Biostrings::readDNAStringSet("protaxFungi/addedmodel/sintaxits2train.fa")
  dbseqs <- dbseqs[grepl(group, names(dbseqs))]
  lengths <-  c(10000, 3000, 1000, 300, 100, 10)
  lengths <- lengths[lengths <= length(dbseqs)]
  # times[[group]] <- tibble::tibble(threads = integer(), minsplit = integer(), length = integer(), time = numeric())
  for (l in lengths) {
    cat(group, l, "\n")
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
    set1 <- sort(sample.int(length(dbseqs), l %/% 2, replace = FALSE))
    dbseqs1 <- dbseqs[set1]
    dbseqs2 <- dbseqs[-set1]
    cat(" clustering with parallel multi-tree method...")
    slmp_time <- system.time(
       slmp <- single_linkage_multi(
          hitfile,
          names(dbseqs),
          dmin = 0.000,
          dmax = 0.4,
          dstep = 0.001,
          preclust = list(names(dbseqs), names(dbseqs1), names(dbseqs2)),
          threads = 8
       )
    )
    cat(" finished in", slmp_time[3], "s\n")
    cat(" clustering with serial multi-tree method...")
    slms_time <- system.time(
       slms <- single_linkage_multi(
          hitfile,
          names(dbseqs),
          dmin = 0.000,
          dmax = 0.4,
          dstep = 0.001,
          preclust = list(names(dbseqs), names(dbseqs1), names(dbseqs2)),
          threads = 1
       )
    )
    cat(" finished in", slms_time[3], "s\n")
    hitfile1 <- sprintf("temp/%s_ref_%d_1.hits", group, l)
    hit_conn1 <- file(hitfile1, open = "w")
    hitfile2 <- sprintf("temp/%s_ref_%d_2.hits", group, l)
    hit_conn2 <- file(hitfile2, open = "w")

    hit_conn = file(hitfile, open = "r")
    dm <- as.dist(matrix(rep(1,l^2), nrow = l))
    cat(" generating distance matrix...")
    n <- 0
    n_write <- 0
    dm_time <- system.time(
       while(length(line <- readLines(hit_conn, 1))) {
          n = n + 1
          distline <- as.numeric(strsplit(line, split = "\t")[[1]])
          if (distline[1] == distline[2]) next
          if (distline[1] < distline[2]) {
             i = distline[1] + 1
             j = distline[2] + 1
          } else {
             j = distline[1] + 1
             i = distline[2] + 1
          }
          il = intersect_length(c(i,j), set1)
          if (il == 2L) {
             n_write = n_write + 1
             writeLines(paste(sum(set1 < i), sum(set1 < j), distline[3], sep = "\t"), hit_conn1)
          } else if (il == 0L) {
             n_write = n_write + 1
             writeLines(paste(i - sum(set1 < i) - 1, j - sum(set1 < j) - 1, distline[3], sep = "\t"), hit_conn2)
          }
          dm[l*(i-1) - i*(i-1)/2 + j-i] <- distline[3]
       }
    )
    cat(" read", n, "distances and wrote", n_write, "distances in", dm_time[3], "s\n")
    close(hit_conn)
    close(hit_conn1)
    close(hit_conn2)
    cat(" clustering with hclust...")
    hc_time <- system.time(
       hc <- hclust(dm, method = "single") %>%
          lapply(0:400 / 1000, cutree, k = NULL, tree = .) %>%
          lapply(function(x) {
             for (i in seq_along(x)) {if (x[i] < i) x[x >= i] = x[x >= i] + 1L}
             x
          }) %>%
          do.call(rbind, .)
    )
    hc <- hc - 1L
    colnames(hc) <- names(dbseqs)
    cat(" finished in", hc_time[3], "s")
    hc_time1 <- system.time(
       hc1 <- hclust(as.dist(as.matrix(dm)[set1,set1]), method = "single") %>%
          lapply(0:400 / 1000, cutree, k = NULL, tree = .) %>%
          lapply(function(x) {
             for (i in seq_along(x)) {if (x[i] < i) x[x >= i] = x[x >= i] + 1L}
             x
          }) %>%
          do.call(rbind, .)
    )
    hc1 <- hc1 - 1L
    colnames(hc1) <- names(dbseqs1)
    cat(" +", hc_time1[3], "s")
    hc_time2 <- system.time(
       hc2 <- hclust(as.dist(as.matrix(dm)[-set1,-set1]), method = "single") %>%
          lapply(0:400 / 1000, cutree, k = NULL, tree = .) %>%
          lapply(function(x) {
             for (i in seq_along(x)) {if (x[i] < i) x[x >= i] = x[x >= i] + 1L}
             x
          }) %>%
          do.call(rbind, .)
    )
    hc2 <- hc2 - 1L
    colnames(hc2) <- names(dbseqs2)
    cat(" +", hc_time2[3], "s =", hc_time[3] + hc_time1[3] + hc_time2[3], "s\n")
    cat(" clustering with matrix method...")
    slmt_time <- system.time(
       slmt <- single_linkage_matrix_thread(
          hitfile,
          names(dbseqs),
          dmin = 0.000,
          dmax = 0.4,
          dstep = 0.001,
          threads = 8
       )
    )
    colnames(slmt) <- names(dbseqs)
    cat(" finished in", slmt_time[3], "s")
    slmt_time1 <- system.time(
       slmt1 <- single_linkage_matrix_thread(
          hitfile1,
          names(dbseqs1),
          dmin = 0.000,
          dmax = 0.4,
          dstep = 0.001,
          threads = 8
       )
    )
    colnames(slmt1) <- names(dbseqs1)
    cat(" +", slmt_time1[3], "s")
    slmt_time2 <- system.time(
       slmt2 <- single_linkage_matrix_thread(
          hitfile2,
          names(dbseqs2),
          dmin = 0.000,
          dmax = 0.4,
          dstep = 0.001,
          threads = 2
       )
    )
    colnames(slmt2) <- names(dbseqs2)
    cat(" +", slmt_time2[3], "s =", slmt_time[3] + slmt_time1[3] + slmt_time2[3], "s\n")
    cat(" clustering with tree method...")
    slp_time <- system.time(
       slp <- single_linkage_pool(
          hitfile,
          names(dbseqs),
          dmin = 0.000,
          dmax = 0.4,
          dstep = 0.001
       )
    )
    colnames(slp) <- names(dbseqs)
    cat(" finished in", slp_time[3], "s")
    slp_time1 <- system.time(
       slp1 <- single_linkage_pool(
          hitfile1,
          names(dbseqs1),
          dmin = 0.000,
          dmax = 0.4,
          dstep = 0.001
       )
    )
    colnames(slp1) <- names(dbseqs1)
    cat(" +", slp_time1[3], "s")
    slp_time2 <- system.time(
       slp2 <- single_linkage_pool(
          hitfile2,
          names(dbseqs2),
          dmin = 0.000,
          dmax = 0.4,
          dstep = 0.001
       )
    )
    colnames(slp2) <- names(dbseqs2)
    cat(" +", slp_time2[3], "s =", slp_time[3] + slp_time1[3] + slp_time2[3], "s\n")
    stopifnot(identical(slmp, list(hc, hc1, hc2)))
    stopifnot(identical(slmp, slms))
    stopifnot(identical(hc, slmt))
    stopifnot(identical(hc1, slmt1))
    stopifnot(identical(hc2, slmt2))
    stopifnot(identical(hc, slp))
    stopifnot(identical(hc1, slp1))
    stopifnot(identical(hc2, slp2))
    cat(" all results match\n")
  }
}
