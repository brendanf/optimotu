
<!-- README.md is generated from README.Rmd. Please edit that file -->

# OptimOTU

<!-- badges: start -->

[![R-CMD-check](https://github.com/brendanf/optimotu/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/brendanf/optimotu/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/brendanf/optimotu/branch/master/graph/badge.svg)](https://app.codecov.io/gh/brendanf/optimotu?branch=master)
<!-- badges: end -->

`optimotu` is a package for cluster optimization, in particular
optimization of thresholds for single-linkage clustering of nucleotide
sequences in order to generate operational taxonomic units (OTUs). It
includes algorithms for single linkage clustering at multiple thresholds
simultaneously, based on a (sparse) distance matrix. Importantly, the
distance matrix is read only once, and not stored in memory, and so it
can be streamed via a linux named pipe, allowing very large sets of
sequences to be efficiently clustered.

Functionality for pairwise distance calculations is in development,
using some of the fastest published libraries for pairwise global
sequence alignment.

## Installation

You can install the development version of optimotu like so:

``` r
install.packages("remotes") # skip if it's already installed
remotes::install_github("brendanf/optimotu")
```

## Example

``` r
library(optimotu)
```

First we download some reference sequences for our taxon of interest and
save them to a file. As an example we will use ITS sequences from
cultures of the genus *Mortierella* which are identified to species
level, generated at the Westerdijk Institute<sup>1</sup>. In real use
cases, assembling and annotating the reference sequences might be the
biggest part of the work!

``` r
morti_search <- rentrez::entrez_search(
  "nuccore",
  term = 'Mortierella[Organism] NOT "Mortierella sp." AND 422523[BioProject] AND internal transcribed spacer',
  retmax = 100
  )
seqs <- rentrez::entrez_fetch("nuccore", id = morti_search$ids, rettype = "fasta")
tf <- tempfile(fileext = ".fasta")
writeLines(seqs, tf)
```

Now we can test cluster the sequences at different thresholds between 0
and 20%, with a step size of 0.1%. Note that for this to work you must
download and install USEARCH!

``` r
clustering <- seq_cluster_usearch(
  seq = tf,
  threshold_config = threshold_uniform(0, 0.2, 0.001, thresh_names = seq(0, 0.2, 0.001))
)
```

The result is a matrix where each column is one of the input sequences,
each row is a threshold, and the entries are the cluster assignment of
the sequence when that clustering threshold is used.

``` r
dim(clustering)
#> [1] 201  55
```

As a convention, `optimotu` always uses the 0-based index of the first
sequence which appears in the cluster as the cluster ID; so the first
sequence is always assigned to cluster 0, the second cluster is assigned
to cluster 0 if it is joined to the first sequence, otherwise it is part
of a new cluster 1, etc.

To assess which threshold is best for dividing the genus *Mortierella*
into species, we need the actual species assignment for each of the
sequences. In this case, the genus and species are the second and third
words in the headers of our downloaded sequences.

``` r
# find the fasta header lines
seq_header <- grep("^>", readLines(tf), value = TRUE)
# split each line into words, and take the second and third words
species <- vapply(
  strsplit(seq_header, " ", fixed = TRUE),
  function(x) paste(x[2], x[3]),
  ""
)

head(species)
#> [1] "Mortierella antarctica"        "Mortierella microzygospora"   
#> [3] "Mortierella alpina"            "Mortierella elongatula"       
#> [5] "Mortierella macrocystopsis"    "Mortierella histoplasmatoides"
```

Now we can calculate scores for how well clustering at different
thresholds matches the true species identities. `optimotu` provides
functions for three types of clustering comparison scores.

The first score type is based in information theory. It consists of the
mutual information and adjusted mutual information<sup>2</sup>. The
`adjusted_mutual_information()` function calculates both of these
metrics, and returns them in a data frame, with the threshold names
passed to `seq_cluster_usearch()` as row names.

``` r
cluster_metrics <- adjusted_mutual_information(clustering, factor(species))
cluster_metrics$threshold <- as.numeric(rownames(cluster_metrics))
```

The second type of score is calculated using a pairwise confusion
matrix, and includes the Rand index, adjusted Rand index, Matthews
correlation coefficient, and Fowlkes-Mallow index<sup>3</sup>. The
confusion matrix approach treats clustering as a binary classification
problem on pairs of sequences, where the classification is considered
“positive” if the two sequences are in the same cluster, and “negative”
if they are in different clusters. These functions can be called on the
test and true clustering (e.g.,
`rand_index(clustering, factor(species))`) but when multiple metrics are
needed, it is more efficient to calculate the confusion matrix once, and
pass it to the functions for each score.

``` r
cm <- confusion_matrix(clustering, factor(species))
cluster_metrics$RI <- rand_index(cm)
cluster_metrics$ARI <- adjusted_rand_index(cm)
cluster_metrics$MCC <- matthews_correlation_coefficient(cm)
cluster_metrics$FMI <- fowlkes_mallow_index(cm)
```

The third type of score is based on cluster-matching. This approach
determines which cluster in the test clustering is the closest match for
each cluster in the true clustering, and scores the test clustering
based on how similar these clusters are. It is the only type of score
which is asymmetric, i.e., the same result is not obtained if the “test”
and “true” clustering are switched. Only the (weighted) F1 score, as
used by Vu et al.<sup>1,4</sup> is included in `optimotu`.

``` r
cluster_metrics$F1 <- fmeasure(clustering, factor(species), 1)
```

Now we can look at how the different scores vary across our range of
clustering thresholds.

``` r
library(ggplot2)
# reshape for easier plotting
cluster_metrics <- 
  tidyr::pivot_longer(
    cluster_metrics,
    cols = c(MI, AMI, RI, ARI, MCC, FMI, F1),
    names_to = "metric",
    values_to = "value"
  )

ggplot(cluster_metrics, aes(x = threshold, y = value, color = metric)) +
  geom_line() +
  facet_wrap(~metric, scales = "free") +
  scale_color_discrete(guide = NULL)
```

<img src="man/figures/README-unnamed-chunk-7-1.png" width="100%" />

Notice that the mutual information (MI) and Rand index (RI) give
comparatively higher scores for thresholds close to 0 (i.e., many small
clusters). This is because these two scores are not adjusted for chance.
Additionally, mutual information is not bounded between 0 and 1. These
two scores are included only for completeness, and it is not recommended
to use them to determine optimum clustering thresholds. The “adjusted”
scores (AMI and ARI, respectively) show clearer peaks away from 0, and
AMI is bounded at 1.

For this case, the confusion-matrix based scores (ARI, MCC, FMI) all
find two peaks of approximately equal magnitude, at 1.9% and 3.7%
sequence dissimilarity, whereas the information theoretic AMI and
cluster-matching F1 both find the peak at 1.9% to be better.

## References

<sup>1</sup> Vu, D., Groenewald, M., de Vries, M., Gehrmann, T.,
Stielow, B., Eberhardt, U., Al-Hatmi, A., Groenewald, J.Z., Cardinali,
G., Houbraken, J., Boekhout, T., Crous, P.W., Robert, V., Verkley,
G.J.M., 2019. Large-scale generation and analysis of filamentous fungal
DNA barcodes boosts coverage for kingdom fungi and reveals thresholds
for fungal species and higher taxon delimitation. Studies in Mycology
92, 135–154. <https://doi.org/10.1016/j.simyco.2018.05.001>

<sup>2</sup> Vinh, N.X., Epps, J., Bailey, J., 2010. Information
Theoretic Measures for Clusterings Comparison: Variants, Properties,
Normalization and Correction for Chance. Journal of Machine Learning
Research 11, 2837–2854.

<sup>3</sup> Warrens, M.J., van der Hoef, H., 2022. Understanding the
Adjusted Rand Index and Other Partition Comparison Indices Based on
Counting Object Pairs. J Classif 39, 487–509.
<https://doi.org/10.1007/s00357-022-09413-z>

<sup>4</sup> Vu, D., Nilsson, R.H., Verkley, G.J.M., 2022. Dnabarcoder:
An open-source software package for analysing and predicting DNA
sequence similarity cutoffs for fungal sequence identification.
Molecular Ecology Resources 22:7.
<https://doi.org/10.1111/1755-0998.13651>
