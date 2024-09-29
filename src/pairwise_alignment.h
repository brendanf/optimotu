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

double distance_wfa2(const std::string&a, const std::string&b, wfa::WFAlignerEdit &aligner);

double distance_edlib(const std::string &a, const std::string &b, EdlibAlignConfig &aligner);

// [[Rcpp::export]]
std::string cigar_wfa2(const std::string &a, const std::string &b,
                       int match = 0, int mismatch = 1,
                       int open1 = 0, int extend1 = 1,
                       int open2 = 0, int extend2 = 1);

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

//' Align two sequences using WFA2 and return the pairwise distance
//'
//' @description WFA2 allows several alignment strategies, and this function
//' selects the least parameterized version possible given the inputs:
//'
//'  - Edit : `gap` == `mismatch` != 0; `match` == `extend` == `gap2` == `extend2` == 0.
//'  - Indel : `gap` != 0; `mismatch` == `match` == `extend` == `gap2` == `extend2` == 0.
//'  - GapLinear : `extend` == `gap2` == `extend2` == 0; other parameters do not meet requirements for "Edit" or "Indel".
//'  - GapAffine : `gap2` == `extend2` == 0; other parameters do not meet requirements for "Edit", "Indel", or "GapLinear".
//'  - GapAffine2Pieces : parameters do not meet requirements for `Edit`, `Indel`, `GapLinear`, or `GapAffine`.
//'
//' @param a (`character` string) first string to align
//' @param b (`character` string) second string to align
//' @param match (`integer` scalar) match score; positive is a bonus.
//' @param mismatch (`integer` scalar) mismatch score; positive is a penalty.
//' @param gap (`integer` scalar) gap opening score; positive is a penalty.
//' Alternatively, if `extend`, `gap2`, and `extend2` are all 0, then this is
//' the penalty for all gaps.
//' @param extend (`integer` scalar) gap extension score; positive is a penalty
//' @param gap2 (`integer` scalar) alternate gap opening score for two-piece
//' affine gap penalty; positive is penalty.  Ignored if both `gap2` and
//' `extend2` are 0.
//' @param extend2 (`integer` scalar) alternate gap extension score for
//' two-piece affine gap penalty; positive is penalty. Ignored if both `gap2`
//' and `extend2` are 0.
//' @export
//'
// [[Rcpp::export]]
double align(const std::string a, const std::string b,
             int match = 0, int mismatch = 1,
             int gap = 1, int extend = 0,
             int gap2 = 0, int extend2 = 0);
#endif //OPTIMOTU_PAIRWISE_ALIGNMENT_H_INCLUDED
