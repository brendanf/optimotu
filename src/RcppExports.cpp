// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppThread.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// confusion_matrix
Rcpp::DataFrame confusion_matrix(const Rcpp::IntegerMatrix k, const Rcpp::IntegerVector c, const int threads);
RcppExport SEXP _optimotu_confusion_matrix(SEXP kSEXP, SEXP cSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type c(cSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(confusion_matrix(k, c, threads));
    return rcpp_result_gen;
END_RCPP
}
// align
double align(std::string a, std::string b);
RcppExport SEXP _optimotu_align(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::string >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(align(a, b));
    return rcpp_result_gen;
END_RCPP
}
// align2
double align2(const std::string& a, const std::string& b);
RcppExport SEXP _optimotu_align2(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(align2(a, b));
    return rcpp_result_gen;
END_RCPP
}
// distmx
Rcpp::DataFrame distmx(std::vector<std::string> seq, double dist_threshold, uint8_t threads, bool heuristic);
RcppExport SEXP _optimotu_distmx(SEXP seqSEXP, SEXP dist_thresholdSEXP, SEXP threadsSEXP, SEXP heuristicSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< double >::type dist_threshold(dist_thresholdSEXP);
    Rcpp::traits::input_parameter< uint8_t >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< bool >::type heuristic(heuristicSEXP);
    rcpp_result_gen = Rcpp::wrap(distmx(seq, dist_threshold, threads, heuristic));
    return rcpp_result_gen;
END_RCPP
}
// distmx_cluster_single
Rcpp::RObject distmx_cluster_single(const std::string file, const Rcpp::CharacterVector seqnames, const Rcpp::List threshold_config, const Rcpp::List method_config, const Rcpp::List parallel_config, const std::string output_type);
RcppExport SEXP _optimotu_distmx_cluster_single(SEXP fileSEXP, SEXP seqnamesSEXP, SEXP threshold_configSEXP, SEXP method_configSEXP, SEXP parallel_configSEXP, SEXP output_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector >::type seqnames(seqnamesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type threshold_config(threshold_configSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type method_config(method_configSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type parallel_config(parallel_configSEXP);
    Rcpp::traits::input_parameter< const std::string >::type output_type(output_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(distmx_cluster_single(file, seqnames, threshold_config, method_config, parallel_config, output_type));
    return rcpp_result_gen;
END_RCPP
}
// distmx_cluster_multi
Rcpp::RObject distmx_cluster_multi(const std::string file, const Rcpp::CharacterVector seqnames, const Rcpp::ListOf<Rcpp::CharacterVector> which, const Rcpp::List threshold_config, const Rcpp::List method_config, const Rcpp::List parallel_config, const std::string output_type);
RcppExport SEXP _optimotu_distmx_cluster_multi(SEXP fileSEXP, SEXP seqnamesSEXP, SEXP whichSEXP, SEXP threshold_configSEXP, SEXP method_configSEXP, SEXP parallel_configSEXP, SEXP output_typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector >::type seqnames(seqnamesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::ListOf<Rcpp::CharacterVector> >::type which(whichSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type threshold_config(threshold_configSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type method_config(method_configSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type parallel_config(parallel_configSEXP);
    Rcpp::traits::input_parameter< const std::string >::type output_type(output_typeSEXP);
    rcpp_result_gen = Rcpp::wrap(distmx_cluster_multi(file, seqnames, which, threshold_config, method_config, parallel_config, output_type));
    return rcpp_result_gen;
END_RCPP
}
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
// fmeasure_list
Rcpp::NumericVector fmeasure_list(Rcpp::ListOf<Rcpp::ListOf<Rcpp::IntegerVector>> k, Rcpp::ListOf<Rcpp::IntegerVector> c, size_t ncpu);
RcppExport SEXP _optimotu_fmeasure_list(SEXP kSEXP, SEXP cSEXP, SEXP ncpuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::ListOf<Rcpp::ListOf<Rcpp::IntegerVector>> >::type k(kSEXP);
    Rcpp::traits::input_parameter< Rcpp::ListOf<Rcpp::IntegerVector> >::type c(cSEXP);
    Rcpp::traits::input_parameter< size_t >::type ncpu(ncpuSEXP);
    rcpp_result_gen = Rcpp::wrap(fmeasure_list(k, c, ncpu));
    return rcpp_result_gen;
END_RCPP
}
// fmeasure_matrix
Rcpp::NumericVector fmeasure_matrix(Rcpp::IntegerMatrix k, Rcpp::IntegerVector c, size_t ncpu);
RcppExport SEXP _optimotu_fmeasure_matrix(SEXP kSEXP, SEXP cSEXP, SEXP ncpuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type c(cSEXP);
    Rcpp::traits::input_parameter< size_t >::type ncpu(ncpuSEXP);
    rcpp_result_gen = Rcpp::wrap(fmeasure_matrix(k, c, ncpu));
    return rcpp_result_gen;
END_RCPP
}
// mutual_information
Rcpp::NumericVector mutual_information(const Rcpp::IntegerMatrix k, const Rcpp::IntegerVector c, int threads);
RcppExport SEXP _optimotu_mutual_information(SEXP kSEXP, SEXP cSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(mutual_information(k, c, threads));
    return rcpp_result_gen;
END_RCPP
}
// adjusted_mutual_information
Rcpp::DataFrame adjusted_mutual_information(const Rcpp::IntegerMatrix k, const Rcpp::IntegerVector c, int threads);
RcppExport SEXP _optimotu_adjusted_mutual_information(SEXP kSEXP, SEXP cSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector >::type c(cSEXP);
    Rcpp::traits::input_parameter< int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(adjusted_mutual_information(k, c, threads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_optimotu_confusion_matrix", (DL_FUNC) &_optimotu_confusion_matrix, 3},
    {"_optimotu_align", (DL_FUNC) &_optimotu_align, 2},
    {"_optimotu_align2", (DL_FUNC) &_optimotu_align2, 2},
    {"_optimotu_distmx", (DL_FUNC) &_optimotu_distmx, 4},
    {"_optimotu_distmx_cluster_single", (DL_FUNC) &_optimotu_distmx_cluster_single, 6},
    {"_optimotu_distmx_cluster_multi", (DL_FUNC) &_optimotu_distmx_cluster_multi, 7},
    {"_optimotu_intersect_length", (DL_FUNC) &_optimotu_intersect_length, 2},
    {"_optimotu_intersect_length_string", (DL_FUNC) &_optimotu_intersect_length_string, 2},
    {"_optimotu_inner_fmeasure", (DL_FUNC) &_optimotu_inner_fmeasure, 3},
    {"_optimotu_fmeasure_list", (DL_FUNC) &_optimotu_fmeasure_list, 3},
    {"_optimotu_fmeasure_matrix", (DL_FUNC) &_optimotu_fmeasure_matrix, 3},
    {"_optimotu_mutual_information", (DL_FUNC) &_optimotu_mutual_information, 3},
    {"_optimotu_adjusted_mutual_information", (DL_FUNC) &_optimotu_adjusted_mutual_information, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_optimotu(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
