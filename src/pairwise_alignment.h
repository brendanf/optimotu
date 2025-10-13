#ifndef OPTIMOTU_PAIRWISE_ALIGNMENT_H_INCLUDED
#define OPTIMOTU_PAIRWISE_ALIGNMENT_H_INCLUDED

#include <vector>
#include <string>
#include <utility>
#include <bindings/cpp/WFAligner.hpp>
#include <edlib.h>
#include "alignment_enums.h"

#ifdef OPTIMOTU_R
#include <Rcpp.h>
#endif //OPTIMOTU_R

template<enum AlignmentSpan span = AlignmentSpan::GLOBAL>
double distance_wfa2(const std::string &a, const std::string &b, wfa::WFAligner &aligner);

double distance_edlib(const std::string &a, const std::string &b, EdlibAlignConfig &aligner);

template<enum AlignmentSpan span = AlignmentSpan::GLOBAL>
double distance_edlib(const std::string &a, const std::string &b);

template<enum AlignmentSpan span = AlignmentSpan::GLOBAL>
std::string cigar_wfa2(const std::string &a, const std::string &b,
                       wfa::WFAligner &aligner);


//' @describeIn pairwise_alignment Compute pairwise global alignment CIGAR with WFA2
//' @export
//' @keywords internal
// [[Rcpp::export]]
std::string cigar_wfa2_global(const std::string &a, const std::string &b,
                       int match = 0, int mismatch = 1,
                       int gap_open = 0, int gap_extend = 1,
                       int gap_open2 = 0, int gap_extend2 = 1);

//' @describeIn pairwise_alignment Compute pairwise extension alignment CIGAR with WFA2
//' @export
//' @keywords internal
// [[Rcpp::export]]
std::string cigar_wfa2_extend(const std::string &a, const std::string &b,
                       int match = 0, int mismatch = 1,
                       int gap_open = 0, int gap_extend = 1,
                       int gap_open2 = 0, int gap_extend2 = 1);

template<enum AlignmentSpan span = AlignmentSpan::GLOBAL>
std::pair<double, std::string> distance_and_cigar_wfa2(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner
);

template<enum AlignmentSpan span = AlignmentSpan::GLOBAL>
std::string cigar_edlib(const std::string &a, const std::string &b);

//' @describeIn pairwise_alignment Compute pairwise global alignment CIGAR with Edlib
//' @export
//' @keywords internal
// [[Rcpp::export]]
std::string cigar_edlib_global(const std::string &a, const std::string &b);

//' @describeIn pairwise_alignment Compute pairwise extension alignment CIGAR with Edlib
//' @export
//' @keywords internal
// [[Rcpp::export]]
std::string cigar_edlib_extend(const std::string &a, const std::string &b);

template<enum AlignmentSpan span = AlignmentSpan::GLOBAL>
std::pair<double, std::string> distance_and_cigar_edlib(
    const std::string &a,
    const std::string &b,
    EdlibAlignConfig &aligner
);

template<enum AlignmentSpan span = AlignmentSpan::GLOBAL>
std::pair<int, double> score_and_distance_wfa2(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner
);

//' @describeIn pairwise_alignment Compute pairwise alignment distance with WFA2
//' @export
//' @keywords internal
// [[Rcpp::export]]
double align_wfa2_global(const std::string a, const std::string b,
             int match = 0, int mismatch = 1,
             int gap_open = 0, int gap_extend = 1,
             int gap_open2 = 0, int gap_extend2 = 0);

//' @describeIn pairwise_alignment Compute pairwise alignment distance with WFA2
//' @export
//' @keywords internal
// [[Rcpp::export]]
double align_wfa2_extend(const std::string a, const std::string b,
                          int match = 0, int mismatch = 1,
                          int gap_open = 0, int gap_extend = 1,
                          int gap_open2 = 0, int gap_extend2 = 0);

//' @describeIn pairwise_alignment Compute pairwise  global alignment distance with Edlib
//' @export
//' @keywords internal
// [[Rcpp::export]]
double align_edlib_global(const std::string a, const std::string b);

//' @describeIn pairwise_alignment Compute pairwise extension alignment distance with Edlib
//' @export
//' @keywords internal
// [[Rcpp::export]]
double align_edlib_extend(const std::string a, const std::string b);

double distance_from_cigar(const std::string &cigar);
double distance_from_cigar_extend(const std::string &cigar);
std::tuple<double, int, int, int, int, int> dist_gapstats_from_cigar(const std::string & cigar);
std::tuple<double, int, int, int, int, int> dist_gapstats_from_cigar_extend(const std::string & cigar);

#ifdef OPTIMOTU_R
//' Add gap statistics to a DataFrame containing CIGAR strings
//' @param df (DataFrame) containing a column with CIGAR strings
//' @param cigar_column (character) name of the column containing the CIGAR strings
//' @return (DataFrame) with the added gap statistics
//' @keywords internal
// [[Rcpp::export]]
Rcpp::RObject add_gapstats(Rcpp::DataFrame df, std::string cigar_column);
#endif //OPTIMOTU_R
#endif //OPTIMOTU_PAIRWISE_ALIGNMENT_H_INCLUDED
