# SPDX-CopyrightText: (c) 2026 Brendan Furneaux
# SPDX-License-Identifier: MIT

# Fraction of estimated swap-free RAM used by clustering_memory_budget_mb =
# "auto". Liberal enough to leave most of the headroom for clustering, while
# still reserving some for R objects, the kernel, and other processes.
AUTO_CLUSTERING_MEMORY_BUDGET_FRACTION <- 0.8

parse_meminfo_kb <- function(key, lines) {
  pattern <- paste0("^", key, ":")
  line <- grep(pattern, lines, value = TRUE)
  if (length(line) != 1L) {
    stop(
      "Could not read ",
      key,
      " from /proc/meminfo.",
      call. = FALSE
    )
  }
  value <- suppressWarnings(
    as.numeric(strsplit(trimws(line), "\\s+")[[1]][2])
  )
  if (!is.finite(value) || value < 0) {
    stop(
      "Could not parse ",
      key,
      " from /proc/meminfo.",
      call. = FALSE
    )
  }
  value
}

cgroup_walk_remaining_bytes <- function(
  fs_root,
  rel_path,
  max_name,
  current_name,
  is_unlimited
) {
  remaining <- Inf
  saw_controller <- FALSE
  if (!nzchar(fs_root)) {
    return(list(remaining = remaining, saw_controller = saw_controller))
  }
  rel_path <- sub("^/", "", rel_path)
  current_dir <- if (!nzchar(rel_path)) {
    fs_root
  } else {
    file.path(fs_root, rel_path)
  }
  fs_root_norm <- normalizePath(fs_root, winslash = "/", mustWork = FALSE)
  seen <- character(0)
  repeat {
    current_norm <- normalizePath(
      current_dir,
      winslash = "/",
      mustWork = FALSE
    )
    if (current_norm %in% seen) {
      break
    }
    seen <- c(seen, current_norm)
    max_file <- file.path(current_dir, max_name)
    if (file.exists(max_file)) {
      max_txt <- readLines(max_file, warn = FALSE, n = 1L)
      if (length(max_txt) == 1L) {
        saw_controller <- TRUE
        if (!isTRUE(is_unlimited(max_txt))) {
          max_bytes <- suppressWarnings(as.numeric(max_txt))
          cur_bytes <- 0
          cur_file <- file.path(current_dir, current_name)
          if (file.exists(cur_file)) {
            cur_txt <- readLines(cur_file, warn = FALSE, n = 1L)
            if (length(cur_txt) == 1L) {
              parsed_cur <- suppressWarnings(as.numeric(cur_txt))
              if (is.finite(parsed_cur) && parsed_cur >= 0) {
                cur_bytes <- parsed_cur
              }
            }
          }
          if (is.finite(max_bytes) && max_bytes >= 0) {
            remaining <- min(remaining, max(0, max_bytes - cur_bytes))
          }
        }
      }
    }
    if (current_norm == fs_root_norm) {
      break
    }
    parent <- dirname(current_dir)
    if (identical(parent, current_dir)) {
      break
    }
    current_dir <- parent
  }
  list(remaining = remaining, saw_controller = saw_controller)
}

linux_cgroup_remaining_bytes <- function(
  cgroup_lines,
  cgroup_fs_root = "/sys/fs/cgroup"
) {
  v2 <- grep("^0::", cgroup_lines, value = TRUE)
  if (length(v2) >= 1L) {
    rel <- sub("^0::", "", v2[[1]])
    walked <- cgroup_walk_remaining_bytes(
      fs_root = cgroup_fs_root,
      rel_path = rel,
      max_name = "memory.max",
      current_name = "memory.current",
      is_unlimited = function(x) identical(x, "max")
    )
    if (isTRUE(walked$saw_controller)) {
      return(walked$remaining)
    }
  }
  v1 <- grep(":[^:]*memory[^:]*:", cgroup_lines, value = TRUE)
  v1 <- v1[!grepl("^0::", v1)]
  if (length(v1) >= 1L) {
    parts <- strsplit(v1[[1]], ":", fixed = TRUE)[[1]]
    rel <- if (length(parts) >= 3L) {
      paste(parts[-c(1L, 2L)], collapse = ":")
    } else {
      ""
    }
    walked <- cgroup_walk_remaining_bytes(
      fs_root = file.path(cgroup_fs_root, "memory"),
      rel_path = rel,
      max_name = "memory.limit_in_bytes",
      current_name = "memory.usage_in_bytes",
      is_unlimited = function(x) {
        n <- suppressWarnings(as.numeric(x))
        !is.finite(n) || n >= 2^62
      }
    )
    if (isTRUE(walked$saw_controller)) {
      return(walked$remaining)
    }
  }
  Inf
}

estimate_auto_clustering_memory_budget_mb <- function(
  sysname = Sys.info()[["sysname"]],
  meminfo_lines = NULL,
  cgroup_lines = NULL,
  cgroup_fs_root = "/sys/fs/cgroup",
  fraction = AUTO_CLUSTERING_MEMORY_BUDGET_FRACTION,
  min_mb = 1
) {
  if (!identical(sysname, "Linux")) {
    stop(
      "`clustering_memory_budget_mb = \"auto\"` is only available on Linux.",
      call. = FALSE
    )
  }
  if (is.null(meminfo_lines)) {
    if (!file.exists("/proc/meminfo")) {
      stop(
        "Could not read /proc/meminfo for auto memory budget.",
        call. = FALSE
      )
    }
    meminfo_lines <- readLines("/proc/meminfo", warn = FALSE)
  }
  avail_mb <- parse_meminfo_kb("MemAvailable", meminfo_lines) / 1024
  if (is.null(cgroup_lines)) {
    cgroup_file <- "/proc/self/cgroup"
    cgroup_lines <- if (file.exists(cgroup_file)) {
      readLines(cgroup_file, warn = FALSE)
    } else {
      character()
    }
  }
  cgroup_bytes <- linux_cgroup_remaining_bytes(
    cgroup_lines = cgroup_lines,
    cgroup_fs_root = cgroup_fs_root
  )
  cgroup_mb <- cgroup_bytes / (1024 * 1024)
  usable_mb <- min(avail_mb, cgroup_mb)
  max(min_mb, usable_mb * fraction)
}

assert_clustering_memory_budget_mb <- function(x, lower = 0) {
  if (is.null(x) || identical(x, "auto")) {
    return(invisible(x))
  }
  if (is.character(x)) {
    stop(
      "`clustering_memory_budget_mb` must be NULL, \"auto\", or a finite ",
      "number.",
      call. = FALSE
    )
  }
  checkmate::assert_number(
    x,
    lower = lower,
    finite = TRUE,
    .var.name = "clustering_memory_budget_mb"
  )
  invisible(x)
}

resolve_clustering_memory_budget_mb <- function(
  x,
  lower = 0,
  ...
) {
  assert_clustering_memory_budget_mb(x, lower = lower)
  if (is.null(x)) {
    return(NULL)
  }
  if (identical(x, "auto")) {
    return(estimate_auto_clustering_memory_budget_mb(...))
  }
  x
}
