library(magrittr)
source("scripts/00_util.R")
source("scripts/01_functions.R")
times <- list()
groups <- c("Fungi", "Penicillium")
for (group in groups) {
  dbseqs <- Biostrings::readDNAStringSet("protaxFungi/addedmodel/sintaxits2train.fa")
  dbseqs <- dbseqs[grepl(group, names(dbseqs))]
  lengths <-  c(3000, 1000, 300, 100, 10)
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
    dm <- as.dist(matrix(rep(1,l^2), nrow = l))
    cat(group, l, "...\n generating distance matrix...")
    n <- 0
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
          dm[l*(i-1) - i*(i-1)/2 + j-i] <- distline[3]
       }
    )
    cat(" read", n, "distances in", dm_time[3], "s\n")
    close(hit_conn)
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
    cat(" finished in", hc_time[3], "s\n")
    cat(" clustering with matrix method...")
    slmt_time <- system.time(
       slmt <- single_linkage_matrix_thread(
          hitfile,
          names(dbseqs),
          dmin = 0.000,
          dmax = 0.4,
          dstep = 0.001
       )
    )
    cat(" finished in", slmt_time[3], "s\n")
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
    cat(" finished in", slp_time[3], "s\n")
    stopifnot(identical(slmt, hc))
    stopifnot(identical(slmt, slp))
    cat(" all results match\n")
  }
}
