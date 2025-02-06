# optimotu (development version)

* Compilation in Windows is now working in at least some circumstances (GitHub
Actions tests)
* Update to WFA2-lib v2.3
* New interface to `align()`, which now takes a `method` argument to choose
  WFA2-lib or Edlib for alignment. This may break old code, because the
  naming of the arguments for gap penalties has changed. Also note that linear
  gap penalties are now specified using `gap_extend` (the new version of
  `extend`) instead of `gap_open` (the new version of `gap`).
* `seq_cluster()` now allows integer values for `verbose`, and produces very
  verbose output with `verbose=2`.
* Added top-level user function `optimize_thresholds()` to optimize thresholds
  across ranks for a set of taxonomically identified reference sequences.
  `calculate_cluster_measures()` calculates multiple clustering quality measures
  for pre-calculated clusters and returns them in a tidy `data.frame`, suitable
  for plotting with `ggplot2::ggplot`.
* Added `clean_taxonomy()`, backed by `is_placeholder()` to remove common 
  placeholder names (such as "*_sp", "*_incertae_sedis", "uncultured *", etc.) from
  a taxonomy table.  This is recommended before using `optimize_thresholds()`.
* Fixed a bug when using `seq_cluster()` or `distmx_cluster()` with a
  `which` argument that is a list (i.e. simultaneous clustering of subsets).
* `clust_tree()` now has arguments `verbose` and `test` to control verbosity
  and whether to run verification tests during clustering. Tests slow down
  clustering, especially at `test = 2`.
* changed the gap open penalties in `seq_cluster_usearch()` to open=0 (was 1)
  for consistency with edlib and the defaults for WFA2.
* add `seq_search()`, which supports all the same distance calculation options
  as `seq_distmx()` and `seq_cluster()`, but returns only the best match(es) for
  each query sequence.
* Added `optimotu()`, a top-level function to run the OptimOTU clustering
  phase, starting from taxonomically identified sequences and the thresholds
  produced by `optimize_thresholds()`.

# optimotu 0.8.1

* Fix compilation in Linux/Mac

# optimotu 0.8.0

* Unify the interface of `seq_clust()` and `seq_clust_usearch()`; `seq_clust()`
now accepts the `which` argument, and delegates to `seq_clust_usearch()` via the
new `dist_usearch()` helper.

# optimotu 0.7.0

* Introduce functions to calculate pairwise distances and output a distance
matrix.  These use edlib and WFA2 for pairwise alignment, and optionally do
"prealignment" filters to reduce the number of full alignments that must be
performed, using EdLib, WFA2 in edit-distance mode, SneakySnake, or a Kmer
filter. The functions are `seq_distmx_edlib()` (EdLib, no filter);
`seq_distmx_wfa2()` (WFA2 alignment, optional WFA2 edit distance filter);
`seq_distmx_snsn()` (WFA2 alignment, SneakySnake filter);
`seq_distmx_kmer()` (WFA2 alignment, Kmer distance filter);
`seq_distmx_hybrid()` (SneakySnake filter, which routes to EdLib or WFA2
alignment depending on predicted distance).  In a future version these will
be combined into a single function with options for different methods
* Introduce `seq_clust()` function, which calculates pairwise sequence distances
followed by clustering.  Options for sequence distance calculation include
published libraries edlib and WFA2, as well as fast Hamming distance calculation
based on code in ProtaxA. Filters are not yet implemented.
* Add SLINK as a clustering method for `seq_clust()`, `distmx_clust()`, and
`distmx_clust_usearch()` with `cluster_slink()`.  This classic method is optimal
when the distances can be guaranteed to be calculated in order.


# optimotu 0.6.5

* Fixed incorrect end-of-file behavior in `distmx_cluster()` which caused
occasional errors when clustering in parallel.

# optimotu 0.6.4

* Fixed incorrect bounds checking in `confusion_matrix()` and `fmeasure()` which
led to occasional index-out-of-bounds errors or segmentation faults (and
presumably some incorrect results).

# optimotu 0.6.3

* Fixed bug which caused an error when using named thresholds with subset
clustering.

# optimotu 0.6.2

* Fixed missing images in README.
* Fixed bugs which cause methods "index" and "tree" to fail in rare cases.
* Removed (internal) buffering for parallel cluster algorithms; this improves
concurrency when clustering subsets (with argument "which").

# optimotu 0.6.1

* Added error checking to detect possible infinite loops in "tree" algorithm.

# optimotu 0.6.0

* New backend C++ code unifies methods for various algorithms for clustering
based on a distance matrix.  All of these are now accessed through one
user-visible function, `distmx_cluster()`.
* Added a new "index"ed matrix clustering method, which is the default for
`distmx_cluster()`.
* Some options which were previously only available for the "tree" clustering
method (i.e., hclust output and subset clustering with the "which" option) are
now available for the "matrix" clustering method, as well as the new "index"
method.
* Speed of the "matrix" clustering method has been improved.
* All clustering methods can be run in three multithreaded modes, controlled by
the option "parallel_config". The "concurrent" mode uses read-write locks on the
data structure to ensure that only one thread writes at the same time. The
"merge" mode keeps separate copies of the data structure, which are merged once
the full input has been processed.  The "hierarchical" mode uses several
copies of the data structure, each of which is operated on a small number of
threads in "concurrent" mode; the different copies are then merged as in "merge"
mode. The benefits of parallel processing for a sparse distance matrix stored
in a file or read from a pipe are fairly minor, since only one thread can read
from the input at a time.
* `threshold_config()` and helpers have a "thresh_names" argument.  These names
are used as rownames in clustering output matrices. This behavior was previously
only implemented in `seq_clust_usearch()`.
* Test builds now work on Windows and OSX in addition to Linux.

# optimotu 0.5.1

* Fix bug which prevented installation on systems without gproftools installed.

# optimotu 0.5.0

* Dramatically improve speed of `adjusted_mutual_information()`, while also
 avoiding numerical underflows.
* Add documentation of clustering quality metrics.

# optimotu 0.4.4

* Fix even more bugs which caused `usearch_single_linkage()` to fail when argument
 'which' was a list.

# optimotu 0.4.3

* Fix additional bugs which caused `usearch_single_linkage()` to fail when argument
 'which' was a list.

# optimotu 0.4.2

* Fix a bug which caused `usearch_single_linkage()` to fail when argument
 'which' was a list.
* Fix a bug in `rand_index()`, `adjusted_rand_index()`,
  `matthews_correlation_coefficient()`, `fowlkes_mallow_index()` which caused
  them to fail when called directly on a cluster matrix instead of on the
  output of `confusion_matrix()`.

# optimotu 0.4.1

* Fix typo in new `fmeasure()` implementation.

# optimotu 0.4.0

* New functions to calculate additional clustering quality metrics:
  `mutual_information()`, `adjusted_mutual_information()`, `confusion_matrix()`,
  `rand_index()`, `adjusted_rand_index`, `matthews_correlation_coefficient()`,
  `fowlkes_mallow_index()`.
* `fmeasure()` now accepts the same matrix representation of clusters as the
  other clustering algorithms (and which is the output of `single_linkage()`).
  This representation leads to much faster computation, but the old
  representation is also still supported.

# optimotu 0.3.0

* `single_linkage()` and `usearch_single_linkage()` now allow explicit threshold
  specifications which are non-decreasing, instead of strictly increasing.

# optimotu 0.2.5

* More bug fixes
* `usearch_single_linkage()` now checks for the presence of `usearch`.

# optimotu 0.2.4

* bug fix

# optimotu 0.2.3

* Bug fix

# optimotu 0.2.2

* Bug fix

# optimotu 0.2.1

* Bug fix

# optimotu 0.2.0

* Added a `NEWS.md` file to track changes to the package.
* Fixed bug when some preclusters are empty.
* **BREAKING CHANGE**: former functions `single_linkage_pool()`,
  `single_linkage_multi()` and `single_linkage_matrix_thread()` are now combined
  in `single_linkage()`, with the argument `method="tree"` for the 'pool' and
  'multi' (with `which`) version, and `method="matrix"` for the 'matrix_thread'
  version. The names of some arguments have also changed (e.g. `dmin` ->
  `thresh_min`).
* Instead of arguments `thresh_min`, `thresh_max`, and `thresh_step`,
  `single_linkage()` can now alternatively take the argument `thresholds`, to
  explicitly define the desired thresholds.  This is slightly slower in cases
  where the desired thresholds are evenly spaced, but allows more flexibility
  when they are not.
* Added function `usearch_single_linkage()` which calculates a distance matrix
  between DNA sequences using the external program usearch, and then clusters
  the result with `single_linkage()`.

# optimotu 0.1.0

* First numbered version
