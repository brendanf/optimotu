# optimotu: Agent Onboarding

This document is for agents working inside the `optimotu` package repository.
It focuses on package-level architecture and maintenance conventions.

## 1) What `optimotu` is responsible for

`optimotu` is the algorithmic core for clustering and threshold optimization.
It provides:

- distance computation backends (`wfa2`, `edlib`, `hamming`, `usearch`, etc.)
- single-linkage clustering across many thresholds
- threshold optimization against known taxonomy
- clustering quality metrics (MCC, ARI, AMI, FMI, FM, etc.)

This package is intentionally lower-level than `optimotu.pipeline`.

## 2) Key module layout

- `R/config.R`
  - user-facing config objects for distance, clustering, thresholds, and
    parallel behavior
- `R/seq_distmx.R`
  - distance matrix/search entry points and backend dispatch
- `R/seq_cluster.R`
  - sequence clustering entry points and method dispatch
- `R/optimize_thresholds.R`
  - threshold optimization and measure-based selection
- `R/confusion_matrix.R`, `R/fmeasure.R`, `R/util.R`
  - metric internals and shared utilities

Native code:

- `src/`
  - performance-critical C/C++ implementations
  - distance workers, clustering algorithms, generators, and glue via Rcpp
  - R package builds define `OPTIMOTU_R` (`src/Makevars`). All R / Rcpp / R
    API use must be wrapped in `#ifdef OPTIMOTU_R`. Core algorithms keep a
    pure C++ entry (standard types / `SequenceView`) so a standalone C++
    library or non-R wrappers remain possible. Dual factories in
    `src/config.h` (C++ types first, Rcpp config overloads under
    `OPTIMOTU_R`) are the pattern to copy.
  - `src/SequenceView.h` is an R-free non-owning view of sequence bytes;
    `sequence_views_from_r()` (R-only) aliases R `CHARSXP`s without copying.
  - `src/MemoryBudget.h` provides best-effort accounting for clustering-owned
    allocations (including Hamming packed sequences), plumbed through
    `config.*`, `MultipleClusterAlgorithm.*`, `seq_cluster.cpp`,
    `cluster_run.cpp`, and the cluster workers. Exceeding the budget raises
    `MemoryBudgetExceeded`, which surfaces in R as an error whose message
    begins `clustering memory budget exceeded`.
    `clustering_memory_budget_mb = "auto"` (non-default, Linux only) sets
    that budget to 80% of `MemAvailable`, capped by remaining cgroup memory
    when a limit is set (`R/memory_budget.R`). Other platforms error.
  - `SingleClusterAlgorithm::write_threshold_row()` fills one threshold's
    cluster-assignment vector without an `n x m` matrix. Tree and SLINK use
    this in `optimize_thresholds()` so the result matrix is never allocated.
  - `optimize_thresholds()` reorders sequences into lexical taxonomy order
    (ranks, then ID) before clustering so Hamming/`AllPairGenerator` merge
    tiles see taxa as contiguous index blocks. Permutation is of character
    vector entries only (no CHARSXP copy). Skipped for
    `dist_file(by_name = FALSE)`.
  - `MultipleClusterAlgorithm::make_child()` builds tile/copy children
    without holding the parent MCA mutex for the whole construction; only
    `children.push_back` is locked so `parallel_merge` tile init can run
    concurrently (per-subset `make_child` locks remain).

## 3) Public API surface (high-impact functions)

Core exported API families include:

- config helpers:
  - `dist_config()` + `dist_wfa2()`, `dist_edlib()`, `dist_hamming()`,
    `dist_usearch()`, `dist_file()`, `dist_hybrid()`
  - `clust_config()` + `clust_tree()`, `clust_matrix()`, `clust_index()`,
    `clust_slink()`
  - `threshold_config()` + `threshold_uniform()`, `threshold_set()`,
    `threshold_lookup()`
  - `parallel_config()` + `parallel_concurrent()`, `parallel_merge()`,
    `parallel_hierarchical()`
- execution:
  - `seq_distmx()`, `seq_cluster()`, `closed_ref_cluster()`, `distmx_cluster()`
- optimization and scoring:
  - `optimize_thresholds()`, `calc_taxon_thresholds()`,
    `calc_subtaxon_thresholds()`, `find_best_threshold()`,
    `calculate_cluster_measures()`

When changing behavior, prioritize preserving these interfaces unless the user
explicitly requests API changes.

## 4) Typical edit zones

- tune backend behavior or defaults:
  - `R/config.R`
- adjust clustering or distance dispatch:
  - `R/seq_cluster.R`, `R/seq_distmx.R`
- change optimization logic:
  - `R/optimize_thresholds.R`
- change native algorithm behavior/performance:
  - `src/*.cpp`, `src/*.h`

## 5) Testing workflow

- Unit tests are in `tests/testthat/`.
- Native code changes should be validated with relevant tests touching altered
  paths (distance, clustering, thresholding, metrics).
- `test-seq_cluster_parallel.R` now includes larger (~25 sequence) synthetic
  correctness coverage (including `clust_tree`, `clust_slink`, and parallel
  modes), plus backend distance consistency checks for close pairs.
- `test-seq_cluster_multi.R` checks that clustering several overlapping
  subsets in one `seq_cluster()` call matches clustering each subset on its
  own, including crossing (non-nested) subsets, `parallel_merge()`, and
  `parallel_concurrent()` with multiple threads.
  `test-seq_cluster_parallel.R` builds overlapping subsets in a fixture but
  only calls `seq_cluster(..., which = TRUE)`, so it does not cover this.
- Multi-subset `parallel_merge()` tile children share parent routing via
  `tile_routing_parent`. Subset shards (tree/matrix/index/SLINK) use
  `subset_tile_fwd_map()` intersection children sized to `|subset ∩ tile|`.
  Parent-subset-local ids are the rank of each sequence in global input
  order; `write_to_matrix()` / `as_hclust()` permute back to the original
  `which` order. Tree/matrix and SLINK merge intersection shards are
  covered in `test-seq_cluster_multi.R`, including shuffled `which` and
  crossing subsets with `parallel_merge(2)`.
- `MappedClusterAlgorithm` intersection constructors bind `Surrogate` to
  `this->fwd_map` (not the ctor parameter) so merge forwards stable parent
  indices.
- `ClusterSLINK` defers allocation: leaf instances hold the SLINK pointer
  representation; merge parents cache unordered child MST edges and replay
  them through SLINK in `finalize()` (no `ClusterTree` delegate).
- `estimate_bytes(n, role, n_tiles)` takes `ClusterInstanceRole::Leaf` or
  `Parent`. SLINK leaf instances use `n * 8 * 4` bytes; parent peak is
  `n * 8 * 4 + 16 * n * T` where `T` is the merge tile count (`~sqrt(threads)`)
  or `threads` for hierarchical full-copy merge (`T = 1` if unknown).
  Tree/matrix/index ignore `role` / `n_tiles`.
- `estimate_subset_memory_mb()` mirrors native `estimate_bytes()` (tree/SLINK
  `O(n)`, matrix/index `O(nm)`). Merge with `threads = 1` matches concurrent
  (one leaf instance; SLINK allocates pointer arrays only). With more threads,
  merge is one parent root (merge-cache sized) plus leaf-sized tile shards for
  SLINK. `include_result = TRUE` adds `4nm` for tree/SLINK only;
  `optimize_thresholds()` does not need this for tree/SLINK because it scores
  via `write_threshold_row()` instead of materializing the result matrix.
  Batch planning keeps a running integer sequence union (`seq_index` from
  `summarize_by_rank()`) and adds Hamming packed-sequence bytes for that
  union (not the sum of subset sizes) when `dist_hamming()` is used; it does
  not call `unique(unlist(seq_id))` per candidate step.
- `summarize_by_rank()` reuses interned CHARSXPs (no `std::string` ID copies)
  and returns `seq_index` (1-based row indices in the taxonomy) alongside
  `seq_id`. `optimize_thresholds()` passes those indices into clustering
  batches so remapping avoids `match()` on sequence names.
- `test-optimize_thresholds_stream.R` checks that row-wise cluster IDs match
  `write_to_matrix()`, and that streaming optima match
  `seq_cluster()` + `find_best_threshold()`, including duplicate thresholds
  and shuffled `which` order.
- Interrupt regression coverage for `seq_cluster()` lives in
  `test-seq_cluster_interrupt.R` and is environment-guarded (opt-in) to avoid
  flaky default CI runs.
- Header changes do not trigger recompiles, because `R CMD INSTALL` tracks only
  `.cpp` timestamps. After editing anything in `src/*.h`, remove the package's
  own objects (`rm -f src/*.o src/optimotu.so`, keeping the `WFA2-lib/` and
  `edlib/` submodule objects) before rebuilding, or you can get link errors or
  silently stale code.
- Prefer running tests in the project container/environment used by the wider
  OptimOTU ecosystem, not only ad hoc local setups.

## 6) Documentation and style conventions

- Package docs use roxygen2 markdown style.
- Do not hand-edit `NAMESPACE` or `.Rd`; generate via roxygen2.
- Keep comments concise and technical; use them for intent and assumptions.
- For non-trivial behavior changes, add/adjust tests with the same patch set.

## 7) Relationship to sister repos

- `optimotu_targets` orchestrates project-specific `targets` pipelines.
- `optimotu.pipeline` provides workflow glue and wrappers around this package.
- Changes in `optimotu` can propagate into both repos, especially when
  exported APIs or config semantics change.

## 8) Keep this file current

If changes alter architecture, API expectations, or testing conventions in this
document, update this file in the same task.

Common triggers:

- exported function additions/removals/renames
- config object semantic changes
- native backend behavior changes
- major test strategy changes
- shifts in repo responsibilities across the three-project ecosystem

## 9) Suggested reading order

1. `DESCRIPTION`
2. `R/config.R`
3. `R/seq_cluster.R`
4. `R/seq_distmx.R`
5. `R/optimize_thresholds.R`
6. `src/` files relevant to your task
7. matching `tests/testthat/` files
