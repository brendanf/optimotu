meminfo_sample <- c(
  "MemTotal:       16000000 kB",
  "MemFree:         1000000 kB",
  "MemAvailable:     4000000 kB"
)

write_cgroup_files <- function(dir, files) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  for (name in names(files)) {
    writeLines(as.character(files[[name]]), file.path(dir, name))
  }
}

testthat::test_that("parse_meminfo_kb reads MemAvailable", {
  got <- optimotu:::parse_meminfo_kb("MemAvailable", meminfo_sample)
  testthat::expect_equal(got, 4000000)
})

testthat::test_that("parse_meminfo_kb errors when the key is missing", {
  testthat::expect_error(
    optimotu:::parse_meminfo_kb("MemAvailable", "MemTotal: 1 kB"),
    "Could not read MemAvailable"
  )
})

testthat::test_that("auto budget is Linux-only", {
  testthat::expect_error(
    optimotu:::estimate_auto_clustering_memory_budget_mb(
      sysname = "Darwin",
      meminfo_lines = meminfo_sample,
      cgroup_lines = character()
    ),
    "only available on Linux"
  )
  testthat::expect_error(
    optimotu:::resolve_clustering_memory_budget_mb(
      "auto",
      sysname = "Windows"
    ),
    "only available on Linux"
  )
})

testthat::test_that("resolve_clustering_memory_budget_mb passes numbers and NULL", {
  testthat::expect_null(
    optimotu:::resolve_clustering_memory_budget_mb(NULL)
  )
  testthat::expect_equal(
    optimotu:::resolve_clustering_memory_budget_mb(64),
    64
  )
  testthat::expect_error(
    optimotu:::resolve_clustering_memory_budget_mb("AUTO"),
    'must be NULL, "auto", or a finite number'
  )
  testthat::expect_error(
    optimotu:::resolve_clustering_memory_budget_mb(-1),
    "clustering_memory_budget_mb"
  )
})

testthat::test_that("auto budget uses 80% of MemAvailable without cgroup limit", {
  got <- optimotu:::estimate_auto_clustering_memory_budget_mb(
    sysname = "Linux",
    meminfo_lines = meminfo_sample,
    cgroup_lines = character(),
    cgroup_fs_root = tempfile("no-cgroup-")
  )
  testthat::expect_equal(got, 4000000 / 1024 * 0.8)
})

testthat::test_that("auto budget is capped by cgroup v2 remaining", {
  root <- tempfile("cgroup-v2-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  write_cgroup_files(
    root,
    list(memory.max = "max", memory.current = "0")
  )
  write_cgroup_files(
    file.path(root, "job"),
    list(
      memory.max = as.character(100 * 1024 * 1024),
      memory.current = as.character(20 * 1024 * 1024)
    )
  )
  got <- optimotu:::estimate_auto_clustering_memory_budget_mb(
    sysname = "Linux",
    meminfo_lines = meminfo_sample,
    cgroup_lines = "0::/job",
    cgroup_fs_root = root
  )
  # 80 MiB remaining * 0.8, which is tighter than MemAvailable.
  testthat::expect_equal(got, 80 * 0.8)
})

testthat::test_that("cgroup v2 walk uses a parent limit when the leaf is max", {
  root <- tempfile("cgroup-v2-parent-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  write_cgroup_files(
    root,
    list(
      memory.max = as.character(50 * 1024 * 1024),
      memory.current = as.character(10 * 1024 * 1024)
    )
  )
  write_cgroup_files(
    file.path(root, "leaf"),
    list(memory.max = "max", memory.current = "0")
  )
  got <- optimotu:::linux_cgroup_remaining_bytes(
    cgroup_lines = "0::/leaf",
    cgroup_fs_root = root
  )
  testthat::expect_equal(got, 40 * 1024 * 1024)
})

testthat::test_that("auto budget is capped by cgroup v1 remaining", {
  root <- tempfile("cgroup-v1-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  write_cgroup_files(
    file.path(root, "memory", "user.slice"),
    list(
      memory.limit_in_bytes = as.character(80 * 1024 * 1024),
      memory.usage_in_bytes = as.character(16 * 1024 * 1024)
    )
  )
  got <- optimotu:::estimate_auto_clustering_memory_budget_mb(
    sysname = "Linux",
    meminfo_lines = meminfo_sample,
    cgroup_lines = "6:memory:/user.slice",
    cgroup_fs_root = root
  )
  testthat::expect_equal(got, 64 * 0.8)
})

testthat::test_that("cgroup v1 treats huge limits as unlimited", {
  root <- tempfile("cgroup-v1-unlim-")
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  write_cgroup_files(
    file.path(root, "memory"),
    list(
      memory.limit_in_bytes = "9223372036854771712",
      memory.usage_in_bytes = "1000"
    )
  )
  got <- optimotu:::linux_cgroup_remaining_bytes(
    cgroup_lines = "6:memory:/",
    cgroup_fs_root = root
  )
  testthat::expect_equal(got, Inf)
})

testthat::test_that("live Linux auto budget is positive and within MemAvailable", {
  testthat::skip_if_not(identical(Sys.info()[["sysname"]], "Linux"))
  testthat::skip_if_not(file.exists("/proc/meminfo"))
  got <- optimotu:::estimate_auto_clustering_memory_budget_mb()
  avail_mb <- optimotu:::parse_meminfo_kb(
    "MemAvailable",
    readLines("/proc/meminfo", warn = FALSE)
  ) /
    1024
  testthat::expect_true(is.finite(got) && got >= 1)
  testthat::expect_lte(got, avail_mb * 0.8 + 1e-6)
})
