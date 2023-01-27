# optimotu (development version)

* Added a `NEWS.md` file to track changes to the package.
* Fixed bug when some preclusters are empty.
* **BREAKING CHANGE**: former functions `single_linkage_pool()` and
  `single_linkage_matrix_thread()` are now combined in `single_linkage()`, with
  the argument `method="tree"` for the 'pool' version, and `method="matrix"` for
  the 'matrix_thread' version.
* Instead of arguments `dmin`, `dmax`, and `dstep`, `single_linkage()` can now
  alternatively take the argument `thresholds`, to explicitly define the desired
  thresholds.  This is slightly slower in cases where the desired thresholds are
  evenly spaced, but allows more flexibility.

# optimotu 0.1.0

* First numbered version
