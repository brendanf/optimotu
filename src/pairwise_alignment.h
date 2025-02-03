#ifndef OPTIMOTU_PAIRWISE_ALIGNMENT_H_INCLUDED
#define OPTIMOTU_PAIRWISE_ALIGNMENT_H_INCLUDED

#include <vector>
#include <string>
#include <utility>
#include <bindings/cpp/WFAligner.hpp>
#include <edlib.h>

#ifdef OPTIMOTU_R
#include <Rcpp.h>
#endif //OPTIMOTU_R

double distance_wfa2(const std::string &a, const std::string &b, wfa::WFAligner &aligner);

double distance_edlib(const std::string &a, const std::string &b, EdlibAlignConfig &aligner);

//' @return (`character(1)`) CIGAR string
//' @export
//' @keywords internal
//' @describeIn pairwise_alignment Generate alignment CIGAR with WFA2
// [[Rcpp::export]]
std::string cigar_wfa2(const std::string &a, const std::string &b,
                       int match = 0, int mismatch = 1,
                       int gap_open = 0, int gap_extend = 1,
                       int gap_open2 = 0, int gap_extend2 = 1);

//' @describeIn pairwise_alignment Generate alignment CIGAR with Edlib
//' @export
//' @keywords internal
// [[Rcpp::export]]
std::string cigar_edlib(const std::string &a, const std::string &b);

std::pair<int, double> score_and_distance_wfa2(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner
);

std::pair<int, double> score_and_distance_wfa2(
    const std::string &a,
    const std::string &b,
    wfa::WFAlignerEdit &aligner
);

//' @describeIn pairwise_alignment Compute pairwise alignment distance with WFA2
//' @export
//' @keywords internal
// [[Rcpp::export]]
double align_wfa2(const std::string a, const std::string b,
             int match = 0, int mismatch = 1,
             int gap_open = 0, int gap_extend = 1,
             int gap_open2 = 0, int gap_extend2 = 0);

//' @describeIn pairwise_alignment Compute pairwise alignment distance with Edlib
//' @export
//' @keywords internal
// [[Rcpp::export]]
double align_edlib(const std::string a, const std::string b);
#endif //OPTIMOTU_PAIRWISE_ALIGNMENT_H_INCLUDED
