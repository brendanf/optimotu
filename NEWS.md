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
