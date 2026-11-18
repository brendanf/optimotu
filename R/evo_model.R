# Standard PAML amino-acid order, used when a model file has no order line.
.evo_paml_order <- function() {
  strsplit("ARNDCQEGHILKMFPSTWYV", "", fixed = TRUE)[[1L]]
}

.evo_aa3 <- c(
  Ala = "A",
  Arg = "R",
  Asn = "N",
  Asp = "D",
  Cys = "C",
  Gln = "Q",
  Glu = "E",
  Gly = "G",
  His = "H",
  Ile = "I",
  Leu = "L",
  Lys = "K",
  Met = "M",
  Phe = "F",
  Pro = "P",
  Ser = "S",
  Thr = "T",
  Trp = "W",
  Tyr = "Y",
  Val = "V"
)

.evo_model_files <- c(
  mtart = "mtArt.dat",
  mtrev = "mtREV24.dat",
  mtrev24 = "mtREV24.dat",
  mtmam = "mtmam.dat",
  jtt = "JTT.dat",
  wag = "wag.dat",
  lg = "lg.dat"
)

.evo_token_letter <- function(tok) {
  if (nchar(tok) == 1L && tok %in% .evo_paml_order()) {
    return(tok)
  }
  hit <- .evo_aa3[tok]
  if (is.na(hit)) {
    return("")
  }
  unname(hit)
}

.evo_order_from_lines <- function(lines) {
  for (line in lines) {
    toks <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
    toks <- toks[nzchar(toks)]
    if (length(toks) != 20L) {
      next
    }
    letters <- vapply(toks, .evo_token_letter, character(1L), USE.NAMES = FALSE)
    if (all(nzchar(letters)) && length(unique(letters)) == 20L) {
      return(letters)
    }
  }
  .evo_paml_order()
}

#' Read a PAML amino-acid exchangeability file
#'
#' Lower triangle of exchangeabilities, then 20 equilibrium frequencies.
#' An amino-acid order line is used when present; otherwise the standard
#' PAML order `ARNDCQEGHILKMFPSTWYV` is assumed.
#'
#' @param path (`character`) path to a `.dat` file
#' @return list with symmetric `S` and named `pi`
#' @keywords internal
parse_paml_exchangeability <- function(path) {
  checkmate::assert_file_exists(path, access = "r")
  lines <- readLines(path, warn = FALSE)
  text <- paste(lines, collapse = " ")
  toks <- strsplit(text, "[[:space:]]+")[[1L]]
  toks <- toks[nzchar(toks)]
  nums <- suppressWarnings(as.numeric(toks))
  nums <- nums[is.finite(nums)]
  if (length(nums) < 210L) {
    stop("`path` does not contain 190 exchangeabilities and 20 frequencies")
  }
  rates <- nums[seq_len(190L)]
  pi <- nums[191:210]
  states <- .evo_order_from_lines(lines)
  s <- matrix(0, 20L, 20L, dimnames = list(states, states))
  k <- 1L
  for (i in 2:20) {
    for (j in seq_len(i - 1L)) {
      s[i, j] <- rates[[k]]
      k <- k + 1L
    }
  }
  s <- s + t(s)
  if (any(pi <= 0)) {
    stop("equilibrium frequencies must be positive")
  }
  pi <- pi / sum(pi)
  names(pi) <- states
  list(S = s, pi = pi)
}

.evo_named_model <- function(name) {
  key <- tolower(name)
  file <- .evo_model_files[[key]]
  if (is.null(file)) {
    stop(
      "`model` must be one of ",
      paste(unique(names(.evo_model_files)), collapse = ", "),
      ", a PAML file path, or a square exchangeability matrix"
    )
  }
  path <- system.file("extdata", "evo", file, package = "optimotu")
  if (!nzchar(path)) {
    stop("vendored model file '", file, "' was not found")
  }
  parse_paml_exchangeability(path)
}

.evo_from_matrix <- function(model, pi) {
  checkmate::assert_matrix(model, mode = "numeric", nrows = 20L, ncols = 20L)
  dn <- dimnames(model)
  if (is.null(dn) || is.null(dn[[1L]]) || is.null(dn[[2L]])) {
    stop("`model` must have single-letter amino-acid dimnames")
  }
  if (!identical(dn[[1L]], dn[[2L]])) {
    stop("`model` row and column names must match")
  }
  letters <- toupper(dn[[1L]])
  if (any(nchar(letters) != 1L) || !setequal(letters, .evo_paml_order())) {
    stop("`model` dimnames must be the 20 standard amino acids")
  }
  dimnames(model) <- list(letters, letters)
  off <- abs(model - t(model))
  diag(off) <- 0
  scale <- max(1, max(abs(model), na.rm = TRUE))
  if (max(off, na.rm = TRUE) > 1e-6 * scale) {
    stop("`model` must be a symmetric exchangeability matrix")
  }
  model <- (model + t(model)) / 2
  diag(model) <- 0
  if (is.null(pi)) {
    stop("a literal exchangeability matrix requires `pi`")
  }
  checkmate::assert_numeric(pi, len = 20L, any.missing = FALSE, lower = 0)
  if (is.null(names(pi))) {
    stop("`pi` must be named with single-letter amino acids")
  }
  pi_names <- toupper(names(pi))
  if (!setequal(pi_names, letters)) {
    stop("`pi` names must match the exchangeability matrix")
  }
  pi <- pi[match(letters, pi_names)]
  if (any(pi <= 0)) {
    stop("`pi` must be positive")
  }
  pi <- pi / sum(pi)
  names(pi) <- letters
  list(S = model, pi = pi)
}

.resolve_evo_model <- function(model, pi) {
  if (is.matrix(model) || (is.numeric(model) && !is.null(dim(model)))) {
    return(.evo_from_matrix(model, pi))
  }
  checkmate::assert_string(model)
  key <- tolower(model)
  if (key %in% names(.evo_model_files)) {
    out <- .evo_named_model(model)
  } else if (file.exists(model)) {
    out <- parse_paml_exchangeability(model)
  } else {
    stop(
      "`model` must be a known abbreviation, an existing PAML file, ",
      "or an exchangeability matrix"
    )
  }
  if (!is.null(pi)) {
    out <- .evo_from_matrix(out$S, pi)
  }
  out
}

.evo_rate_matrix <- function(s, pi) {
  q <- s * rep(pi, each = length(pi))
  diag(q) <- 0
  diag(q) <- -rowSums(q)
  rate <- -sum(pi * diag(q))
  if (!is.finite(rate) || rate <= 0) {
    stop("exchangeability matrix has a non-positive mean rate")
  }
  q / rate
}

# Q is reversible, so D Q D^{-1} is symmetric for D = diag(sqrt(pi)).
# A general eigen solver does not reconstruct the identity on every
# BLAS (GitHub' CUDA test image does not). The symmetric factorization
# does, and its inverse is U' D rather than a second solve().
.evo_decompose <- function(q, pi) {
  scale <- sqrt(as.numeric(pi))
  balanced <- (scale * q) / rep(scale, each = length(scale))
  balanced <- (balanced + t(balanced)) / 2
  ev <- eigen(balanced, symmetric = TRUE)
  vectors <- ev$vectors
  list(
    lambda = ev$values,
    evec = vectors / scale,
    evec_inv = t(vectors * scale)
  )
}

.evo_time_grid <- function(t_max, n) {
  n <- as.integer(n)
  pos <- exp(seq(log(t_max / 1000), log(t_max), length.out = n - 1L))
  c(0, pos)
}

.evo_logp_luts <- function(states, pi, lambda, evec, evec_inv, t_grid) {
  codes <- aa_code_from_char(states)
  n_grid <- length(t_grid)
  logp <- numeric(n_grid * 32L * 32L)
  for (g in seq_along(t_grid)) {
    p <- evec %*% diag(exp(lambda * t_grid[[g]]), 20L) %*% evec_inv
    p[!is.finite(p) | p < 1e-300] <- 1e-300
    for (i in seq_len(20L)) {
      for (j in seq_len(20L)) {
        val <- log(pi[[i]]) + log(p[i, j])
        idx <- (codes[[i]] * 32L + codes[[j]]) * n_grid + g
        logp[[idx]] <- val
      }
    }
  }
  logp
}

#' @param model (`character` or square `numeric` matrix) an abbreviation
#' (`mtart`, `mtrev`, `mtmam`, `jtt`, `wag`, `lg`), a PAML `.dat` path, or
#' a symmetric exchangeability matrix with single-letter dimnames.
#' Abbreviations and files include published equilibrium frequencies.
#' @param pi (`NULL` or named `numeric`) equilibrium frequencies. Required
#' for a literal matrix. When set with an abbreviation or file, replaces
#' that model's frequencies.
#' @param t_max (`numeric`) largest divergence, in expected substitutions
#' per site. Pairs whose likelihood is still increasing there, and pairs
#' with fewer than `min_overlap` compared sites, return `t_max`.
#' @param t_grid (`integer`) number of likelihood grid points from 0 to
#' `t_max`. Interior peaks are refined by Newton steps.
#' @param min_overlap (`integer`) minimum number of compared residues.
#' @param ignore_gaps (`logical`) retained for the same call pattern as
#' [dist_aa()]. Gaps and nonstandard residues (`X`, `B`, `Z`, `*`) are
#' missing data and are not scored.
#' @param genetic_code (`NULL` or NCBI translation table number) if set,
#' inputs are aligned nucleotides and are translated before scoring.
#' @param codon_offset (`integer` in 0, 1, 2) first codon base in the
#' aligned nucleotides when `genetic_code` is set.
#'
#' @details
#' The distance is the maximum-likelihood divergence under the reversible
#' model `Q_ij = s_ij * pi_j`, scaled so the mean rate is 1. Vendored
#' tables are mtART (Abascal, Posada, and Zardoya 2007), mtREV24 (Adachi
#' and Hasegawa 1996), mtMam (Yang, Nielsen, and Hasegawa 1998), JTT
#' (Jones, Taylor, and Thornton 1992), WAG (Whelan and Goldman 2001), and
#' LG (Le and Gascuel 2008). When the package was built with CUDA and a
#' device is visible, `seq_distmx()`, `seq_search()`, `seq_cluster()`,
#' and `optimize_thresholds()` fill distance tiles on every selected GPU
#' (`OPTIMOTU_CUDA_DEVICES`, default all). Otherwise the CPU Newton path
#' is used with no message.
#'
#' @export
#' @describeIn dist_config helper function for method `"evo"`
dist_evo <- function(
  model = "mtart",
  pi = NULL,
  t_max = 5,
  t_grid = 32L,
  min_overlap = 0L,
  ignore_gaps = TRUE,
  genetic_code = NULL,
  codon_offset = 0L
) {
  checkmate::assert_number(t_max, lower = 1e-6, finite = TRUE)
  checkmate::assert_integerish(
    t_grid,
    lower = 3L,
    upper = 256L,
    len = 1L,
    any.missing = FALSE
  )
  t_grid <- as.integer(t_grid)
  checkmate::assert_count(min_overlap)
  checkmate::assert_flag(ignore_gaps)
  checkmate::assert_integerish(codon_offset, len = 1L, lower = 0L, upper = 2L)
  if (!is.null(genetic_code)) {
    checkmate::assert_integerish(genetic_code, len = 1L, lower = 1L)
    genetic_code <- as.integer(genetic_code)
  }
  resolved <- .resolve_evo_model(model, pi)
  states <- colnames(resolved$S)
  q <- .evo_rate_matrix(resolved$S, resolved$pi)
  ev <- .evo_decompose(q, resolved$pi)
  lambda <- ev$lambda
  evec <- ev$evec
  evec_inv <- ev$evec_inv
  p0 <- evec %*% evec_inv
  if (max(abs(p0 - diag(20L))) > 1e-8) {
    stop("rate matrix eigen decomposition did not reconstruct the identity")
  }
  grid <- .evo_time_grid(t_max, t_grid)
  logp <- .evo_logp_luts(
    states,
    resolved$pi,
    lambda,
    evec,
    evec_inv,
    grid
  )
  model_name <- if (is.character(model) && length(model) == 1L) {
    model
  } else {
    "matrix"
  }
  structure(
    list(
      method = "evo",
      model = model_name,
      min_overlap = as.integer(min_overlap),
      ignore_gaps = ignore_gaps,
      genetic_code = genetic_code,
      codon_offset = as.integer(codon_offset),
      t_max = t_max,
      t_grid = grid,
      exchangeability = resolved$S,
      pi = resolved$pi,
      states = paste(states, collapse = ""),
      lambda = lambda,
      evec = as.numeric(evec),
      evec_inv = as.numeric(evec_inv),
      logp = logp,
      call = match.call()
    ),
    class = "optimotu_dist_config"
  )
}
