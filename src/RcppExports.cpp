// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// intersect_length
int intersect_length(const std::vector<int>& c, const std::vector<int>& k);
RcppExport SEXP _optimotu_intersect_length(SEXP cSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(intersect_length(c, k));
    return rcpp_result_gen;
END_RCPP
}
// intersect_length_string
int intersect_length_string(const std::vector<std::string>& c, const std::vector<std::string>& k);
RcppExport SEXP _optimotu_intersect_length_string(SEXP cSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type c(cSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::string>& >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(intersect_length_string(c, k));
    return rcpp_result_gen;
END_RCPP
}
// inner_fmeasure
double inner_fmeasure(const std::vector<int>& cj, const Rcpp::ListOf<Rcpp::IntegerVector>& kpartition, const std::vector<int>& nk);
RcppExport SEXP _optimotu_inner_fmeasure(SEXP cjSEXP, SEXP kpartitionSEXP, SEXP nkSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<int>& >::type cj(cjSEXP);
    Rcpp::traits::input_parameter< const Rcpp::ListOf<Rcpp::IntegerVector>& >::type kpartition(kpartitionSEXP);
    Rcpp::traits::input_parameter< const std::vector<int>& >::type nk(nkSEXP);
    rcpp_result_gen = Rcpp::wrap(inner_fmeasure(cj, kpartition, nk));
    return rcpp_result_gen;
END_RCPP
}
// fmeasure
Rcpp::NumericVector fmeasure(Rcpp::ListOf<Rcpp::ListOf<Rcpp::IntegerVector>> k, Rcpp::ListOf<Rcpp::IntegerVector> c, size_t ncpu);
RcppExport SEXP _optimotu_fmeasure(SEXP kSEXP, SEXP cSEXP, SEXP ncpuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::ListOf<Rcpp::ListOf<Rcpp::IntegerVector>> >::type k(kSEXP);
    Rcpp::traits::input_parameter< Rcpp::ListOf<Rcpp::IntegerVector> >::type c(cSEXP);
    Rcpp::traits::input_parameter< size_t >::type ncpu(ncpuSEXP);
    rcpp_result_gen = Rcpp::wrap(fmeasure(k, c, ncpu));
    return rcpp_result_gen;
END_RCPP
}
// single_linkage_matrix_thread
Rcpp::IntegerMatrix single_linkage_matrix_thread(const std::string file, const Rcpp::CharacterVector& seqnames, const float dmin, const float dmax, const float dstep, const int threads, const int minsplit);
RcppExport SEXP _optimotu_single_linkage_matrix_thread(SEXP fileSEXP, SEXP seqnamesSEXP, SEXP dminSEXP, SEXP dmaxSEXP, SEXP dstepSEXP, SEXP threadsSEXP, SEXP minsplitSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type seqnames(seqnamesSEXP);
    Rcpp::traits::input_parameter< const float >::type dmin(dminSEXP);
    Rcpp::traits::input_parameter< const float >::type dmax(dmaxSEXP);
    Rcpp::traits::input_parameter< const float >::type dstep(dstepSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< const int >::type minsplit(minsplitSEXP);
    rcpp_result_gen = Rcpp::wrap(single_linkage_matrix_thread(file, seqnames, dmin, dmax, dstep, threads, minsplit));
    return rcpp_result_gen;
END_RCPP
}
// single_linkage_pool
Rcpp::IntegerMatrix single_linkage_pool(const std::string file, const Rcpp::CharacterVector& seqnames, const float dmin, const float dmax, const float dstep);
RcppExport SEXP _optimotu_single_linkage_pool(SEXP fileSEXP, SEXP seqnamesSEXP, SEXP dminSEXP, SEXP dmaxSEXP, SEXP dstepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type seqnames(seqnamesSEXP);
    Rcpp::traits::input_parameter< const float >::type dmin(dminSEXP);
    Rcpp::traits::input_parameter< const float >::type dmax(dmaxSEXP);
    Rcpp::traits::input_parameter< const float >::type dstep(dstepSEXP);
    rcpp_result_gen = Rcpp::wrap(single_linkage_pool(file, seqnames, dmin, dmax, dstep));
    return rcpp_result_gen;
END_RCPP
}
// single_linkage_multi
Rcpp::List single_linkage_multi(const std::string file, const Rcpp::CharacterVector& seqnames, const float dmin, const float dmax, const float dstep, const Rcpp::ListOf<Rcpp::CharacterVector>& preclust, const size_t threads);
RcppExport SEXP _optimotu_single_linkage_multi(SEXP fileSEXP, SEXP seqnamesSEXP, SEXP dminSEXP, SEXP dmaxSEXP, SEXP dstepSEXP, SEXP preclustSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type seqnames(seqnamesSEXP);
    Rcpp::traits::input_parameter< const float >::type dmin(dminSEXP);
    Rcpp::traits::input_parameter< const float >::type dmax(dmaxSEXP);
    Rcpp::traits::input_parameter< const float >::type dstep(dstepSEXP);
    Rcpp::traits::input_parameter< const Rcpp::ListOf<Rcpp::CharacterVector>& >::type preclust(preclustSEXP);
    Rcpp::traits::input_parameter< const size_t >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(single_linkage_multi(file, seqnames, dmin, dmax, dstep, preclust, threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_optimotu_intersect_length", (DL_FUNC) &_optimotu_intersect_length, 2},
    {"_optimotu_intersect_length_string", (DL_FUNC) &_optimotu_intersect_length_string, 2},
    {"_optimotu_inner_fmeasure", (DL_FUNC) &_optimotu_inner_fmeasure, 3},
    {"_optimotu_fmeasure", (DL_FUNC) &_optimotu_fmeasure, 3},
    {"_optimotu_single_linkage_matrix_thread", (DL_FUNC) &_optimotu_single_linkage_matrix_thread, 7},
    {"_optimotu_single_linkage_pool", (DL_FUNC) &_optimotu_single_linkage_pool, 5},
    {"_optimotu_single_linkage_multi", (DL_FUNC) &_optimotu_single_linkage_multi, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_optimotu(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
