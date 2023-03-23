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
