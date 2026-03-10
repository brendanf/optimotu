#include <Rcpp.h>
#include <unordered_set>

//' Verify that `which` is a valid list of sequence subsets
//' @param which (`Rcpp::List`) The list of sequence subsets.
//' @param seqnames (`Rcpp::CharacterVector`) The names of the sequences.
//' @return (`void`) Throws an error if `which` is not a valid list of sequence subsets.
// [[Rcpp::export]]

void verify_which_impl(
  Rcpp::List which,
  const Rcpp::CharacterVector seqnames
) {
  // R strings are stored as pointers, and typically identical strings are
  // pointers to the same *char in a global string pool.
  // Thus we can use a set of pointers to the strings to check for membership.
  std::unordered_set<const char *> seqnames_set;
  for (const auto & seqname : seqnames) {
    seqnames_set.insert(Rcpp::String(seqname).get_cstring());
  }

  int i = 1;
  for (const auto & subset : which) {
    if (!Rcpp::is<Rcpp::CharacterVector>(subset)) {
      Rcpp::stop("Subset is not a character vector");
    }
    const auto & subset_vec = Rcpp::as<Rcpp::CharacterVector>(subset);
    int j = 1;
    for (const auto & seqname : subset_vec) {
      if (seqname == NA_STRING) {
        Rcpp::stop("'which' contains NA, which is not a valid sequence name");
      }
      if (seqnames_set.count(Rcpp::String(seqname).get_cstring()) == 0) {
        // no hit for the pointer
        Rcpp::stop("Sequence name '%s' at position %d in subset %d not found in seqnames", Rcpp::String(seqname).get_cstring(), j, i);
      }
      j++;
    }
    i++;
  }
}