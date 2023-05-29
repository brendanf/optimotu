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

//' @export
//'
// [[Rcpp::export]]
double align(const std::string a, const std::string b,
             int match = 0, int mismatch = 1,
             int gap = 1, int extend = 0,
             int gap2 = 0, int extend2 = 0);
#endif //OPTIMOTU_PAIRWISE_ALIGNMENT_H_INCLUDED
