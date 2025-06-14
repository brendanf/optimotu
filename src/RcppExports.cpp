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
// confusion_matrix2
Rcpp::DataFrame confusion_matrix2(const Rcpp::IntegerMatrix k, const Rcpp::IntegerMatrix c, const int threads);
RcppExport SEXP _optimotu_confusion_matrix2(SEXP kSEXP, SEXP cSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix >::type k(kSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix >::type c(cSEXP);
    Rcpp::traits::input_parameter< const int >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(confusion_matrix2(k, c, threads));
    return rcpp_result_gen;
END_RCPP
}
// distmx_cluster_single
Rcpp::RObject distmx_cluster_single(const std::string file, Rcpp::CharacterVector seqnames, Rcpp::List threshold_config, Rcpp::List clust_config, Rcpp::List parallel_config, const std::string output_type, const bool verbose, const bool by_name);
RcppExport SEXP _optimotu_distmx_cluster_single(SEXP fileSEXP, SEXP seqnamesSEXP, SEXP threshold_configSEXP, SEXP clust_configSEXP, SEXP parallel_configSEXP, SEXP output_typeSEXP, SEXP verboseSEXP, SEXP by_nameSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type file(fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type seqnames(seqnamesSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type threshold_config(threshold_configSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type clust_config(clust_configSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type parallel_config(parallel_configSEXP);
    Rcpp::traits::input_parameter< const std::string >::type output_type(output_typeSEXP);
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool >::type by_name(by_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(distmx_cluster_single(file, seqnames, threshold_config, clust_config, parallel_config, output_type, verbose, by_name));
    return rcpp_result_gen;
END_RCPP
}
// distmx_cluster_multi
Rcpp::List distmx_cluster_multi(const std::string file, const Rcpp::CharacterVector seqnames, const Rcpp::ListOf<Rcpp::CharacterVector> which, const Rcpp::List threshold_config, const Rcpp::List method_config, const Rcpp::List parallel_config, const std::string output_type, const bool verbose, const bool by_name);
RcppExport SEXP _optimotu_distmx_cluster_multi(SEXP fileSEXP, SEXP seqnamesSEXP, SEXP whichSEXP, SEXP threshold_configSEXP, SEXP method_configSEXP, SEXP parallel_configSEXP, SEXP output_typeSEXP, SEXP verboseSEXP, SEXP by_nameSEXP) {
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
    Rcpp::traits::input_parameter< const bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< const bool >::type by_name(by_nameSEXP);
    rcpp_result_gen = Rcpp::wrap(distmx_cluster_multi(file, seqnames, which, threshold_config, method_config, parallel_config, output_type, verbose, by_name));
    return rcpp_result_gen;
END_RCPP
}
// seq_distmx_edlib
Rcpp::DataFrame seq_distmx_edlib(std::vector<std::string> seq, double dist_threshold, bool constrain, std::uint8_t threads);
RcppExport SEXP _optimotu_seq_distmx_edlib(SEXP seqSEXP, SEXP dist_thresholdSEXP, SEXP constrainSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< double >::type dist_threshold(dist_thresholdSEXP);
    Rcpp::traits::input_parameter< bool >::type constrain(constrainSEXP);
    Rcpp::traits::input_parameter< std::uint8_t >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_distmx_edlib(seq, dist_threshold, constrain, threads));
    return rcpp_result_gen;
END_RCPP
}
// fastq_names
Rcpp::CharacterVector fastq_names(std::string x);
RcppExport SEXP _optimotu_fastq_names(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(fastq_names(x));
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
// seq_distmx_kmer
Rcpp::DataFrame seq_distmx_kmer(std::vector<std::string> seq, double dist_threshold, double udist_threshold, int match, int mismatch, int gap_open, int gap_extend, int gap_open2, int gap_extend2, std::uint8_t threads);
RcppExport SEXP _optimotu_seq_distmx_kmer(SEXP seqSEXP, SEXP dist_thresholdSEXP, SEXP udist_thresholdSEXP, SEXP matchSEXP, SEXP mismatchSEXP, SEXP gap_openSEXP, SEXP gap_extendSEXP, SEXP gap_open2SEXP, SEXP gap_extend2SEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< double >::type dist_threshold(dist_thresholdSEXP);
    Rcpp::traits::input_parameter< double >::type udist_threshold(udist_thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type match(matchSEXP);
    Rcpp::traits::input_parameter< int >::type mismatch(mismatchSEXP);
    Rcpp::traits::input_parameter< int >::type gap_open(gap_openSEXP);
    Rcpp::traits::input_parameter< int >::type gap_extend(gap_extendSEXP);
    Rcpp::traits::input_parameter< int >::type gap_open2(gap_open2SEXP);
    Rcpp::traits::input_parameter< int >::type gap_extend2(gap_extend2SEXP);
    Rcpp::traits::input_parameter< std::uint8_t >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_distmx_kmer(seq, dist_threshold, udist_threshold, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, threads));
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
// pairwise_alignment
std::vector<std::string> pairwise_alignment(std::string a, std::string b, Rcpp::List dist_config, int span);
RcppExport SEXP _optimotu_pairwise_alignment(SEXP aSEXP, SEXP bSEXP, SEXP dist_configSEXP, SEXP spanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type a(aSEXP);
    Rcpp::traits::input_parameter< std::string >::type b(bSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type dist_config(dist_configSEXP);
    Rcpp::traits::input_parameter< int >::type span(spanSEXP);
    rcpp_result_gen = Rcpp::wrap(pairwise_alignment(a, b, dist_config, span));
    return rcpp_result_gen;
END_RCPP
}
// cigar_wfa2_global
std::string cigar_wfa2_global(const std::string& a, const std::string& b, int match, int mismatch, int gap_open, int gap_extend, int gap_open2, int gap_extend2);
RcppExport SEXP _optimotu_cigar_wfa2_global(SEXP aSEXP, SEXP bSEXP, SEXP matchSEXP, SEXP mismatchSEXP, SEXP gap_openSEXP, SEXP gap_extendSEXP, SEXP gap_open2SEXP, SEXP gap_extend2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type match(matchSEXP);
    Rcpp::traits::input_parameter< int >::type mismatch(mismatchSEXP);
    Rcpp::traits::input_parameter< int >::type gap_open(gap_openSEXP);
    Rcpp::traits::input_parameter< int >::type gap_extend(gap_extendSEXP);
    Rcpp::traits::input_parameter< int >::type gap_open2(gap_open2SEXP);
    Rcpp::traits::input_parameter< int >::type gap_extend2(gap_extend2SEXP);
    rcpp_result_gen = Rcpp::wrap(cigar_wfa2_global(a, b, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2));
    return rcpp_result_gen;
END_RCPP
}
// cigar_wfa2_extend
std::string cigar_wfa2_extend(const std::string& a, const std::string& b, int match, int mismatch, int gap_open, int gap_extend, int gap_open2, int gap_extend2);
RcppExport SEXP _optimotu_cigar_wfa2_extend(SEXP aSEXP, SEXP bSEXP, SEXP matchSEXP, SEXP mismatchSEXP, SEXP gap_openSEXP, SEXP gap_extendSEXP, SEXP gap_open2SEXP, SEXP gap_extend2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type match(matchSEXP);
    Rcpp::traits::input_parameter< int >::type mismatch(mismatchSEXP);
    Rcpp::traits::input_parameter< int >::type gap_open(gap_openSEXP);
    Rcpp::traits::input_parameter< int >::type gap_extend(gap_extendSEXP);
    Rcpp::traits::input_parameter< int >::type gap_open2(gap_open2SEXP);
    Rcpp::traits::input_parameter< int >::type gap_extend2(gap_extend2SEXP);
    rcpp_result_gen = Rcpp::wrap(cigar_wfa2_extend(a, b, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2));
    return rcpp_result_gen;
END_RCPP
}
// cigar_edlib_global
std::string cigar_edlib_global(const std::string& a, const std::string& b);
RcppExport SEXP _optimotu_cigar_edlib_global(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(cigar_edlib_global(a, b));
    return rcpp_result_gen;
END_RCPP
}
// cigar_edlib_extend
std::string cigar_edlib_extend(const std::string& a, const std::string& b);
RcppExport SEXP _optimotu_cigar_edlib_extend(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::string& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(cigar_edlib_extend(a, b));
    return rcpp_result_gen;
END_RCPP
}
// align_wfa2_global
double align_wfa2_global(const std::string a, const std::string b, int match, int mismatch, int gap_open, int gap_extend, int gap_open2, int gap_extend2);
RcppExport SEXP _optimotu_align_wfa2_global(SEXP aSEXP, SEXP bSEXP, SEXP matchSEXP, SEXP mismatchSEXP, SEXP gap_openSEXP, SEXP gap_extendSEXP, SEXP gap_open2SEXP, SEXP gap_extend2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::string >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type match(matchSEXP);
    Rcpp::traits::input_parameter< int >::type mismatch(mismatchSEXP);
    Rcpp::traits::input_parameter< int >::type gap_open(gap_openSEXP);
    Rcpp::traits::input_parameter< int >::type gap_extend(gap_extendSEXP);
    Rcpp::traits::input_parameter< int >::type gap_open2(gap_open2SEXP);
    Rcpp::traits::input_parameter< int >::type gap_extend2(gap_extend2SEXP);
    rcpp_result_gen = Rcpp::wrap(align_wfa2_global(a, b, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2));
    return rcpp_result_gen;
END_RCPP
}
// align_wfa2_extend
double align_wfa2_extend(const std::string a, const std::string b, int match, int mismatch, int gap_open, int gap_extend, int gap_open2, int gap_extend2);
RcppExport SEXP _optimotu_align_wfa2_extend(SEXP aSEXP, SEXP bSEXP, SEXP matchSEXP, SEXP mismatchSEXP, SEXP gap_openSEXP, SEXP gap_extendSEXP, SEXP gap_open2SEXP, SEXP gap_extend2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::string >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type match(matchSEXP);
    Rcpp::traits::input_parameter< int >::type mismatch(mismatchSEXP);
    Rcpp::traits::input_parameter< int >::type gap_open(gap_openSEXP);
    Rcpp::traits::input_parameter< int >::type gap_extend(gap_extendSEXP);
    Rcpp::traits::input_parameter< int >::type gap_open2(gap_open2SEXP);
    Rcpp::traits::input_parameter< int >::type gap_extend2(gap_extend2SEXP);
    rcpp_result_gen = Rcpp::wrap(align_wfa2_extend(a, b, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2));
    return rcpp_result_gen;
END_RCPP
}
// align_edlib_global
double align_edlib_global(const std::string a, const std::string b);
RcppExport SEXP _optimotu_align_edlib_global(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::string >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(align_edlib_global(a, b));
    return rcpp_result_gen;
END_RCPP
}
// align_edlib_extend
double align_edlib_extend(const std::string a, const std::string b);
RcppExport SEXP _optimotu_align_edlib_extend(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type a(aSEXP);
    Rcpp::traits::input_parameter< const std::string >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(align_edlib_extend(a, b));
    return rcpp_result_gen;
END_RCPP
}
// seq_distmx_wfa2
Rcpp::DataFrame seq_distmx_wfa2(std::vector<std::string> seq, double dist_threshold, int match, int mismatch, int gap_open, int gap_extend, int gap_open2, int gap_extend2, bool prealign, bool constrain, std::uint8_t threads);
RcppExport SEXP _optimotu_seq_distmx_wfa2(SEXP seqSEXP, SEXP dist_thresholdSEXP, SEXP matchSEXP, SEXP mismatchSEXP, SEXP gap_openSEXP, SEXP gap_extendSEXP, SEXP gap_open2SEXP, SEXP gap_extend2SEXP, SEXP prealignSEXP, SEXP constrainSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< double >::type dist_threshold(dist_thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type match(matchSEXP);
    Rcpp::traits::input_parameter< int >::type mismatch(mismatchSEXP);
    Rcpp::traits::input_parameter< int >::type gap_open(gap_openSEXP);
    Rcpp::traits::input_parameter< int >::type gap_extend(gap_extendSEXP);
    Rcpp::traits::input_parameter< int >::type gap_open2(gap_open2SEXP);
    Rcpp::traits::input_parameter< int >::type gap_extend2(gap_extend2SEXP);
    Rcpp::traits::input_parameter< bool >::type prealign(prealignSEXP);
    Rcpp::traits::input_parameter< bool >::type constrain(constrainSEXP);
    Rcpp::traits::input_parameter< std::uint8_t >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_distmx_wfa2(seq, dist_threshold, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, prealign, constrain, threads));
    return rcpp_result_gen;
END_RCPP
}
// seq_cluster_single
Rcpp::RObject seq_cluster_single(const Rcpp::CharacterVector& seq, const Rcpp::List dist_config, const Rcpp::List threshold_config, const Rcpp::List clust_config, const Rcpp::List parallel_config, const std::string output_type, const int verbose);
RcppExport SEXP _optimotu_seq_cluster_single(SEXP seqSEXP, SEXP dist_configSEXP, SEXP threshold_configSEXP, SEXP clust_configSEXP, SEXP parallel_configSEXP, SEXP output_typeSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type dist_config(dist_configSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type threshold_config(threshold_configSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type clust_config(clust_configSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type parallel_config(parallel_configSEXP);
    Rcpp::traits::input_parameter< const std::string >::type output_type(output_typeSEXP);
    Rcpp::traits::input_parameter< const int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_cluster_single(seq, dist_config, threshold_config, clust_config, parallel_config, output_type, verbose));
    return rcpp_result_gen;
END_RCPP
}
// seq_cluster_multi
Rcpp::List seq_cluster_multi(const Rcpp::CharacterVector& seq, const Rcpp::ListOf<Rcpp::CharacterVector> which, const Rcpp::List dist_config, const Rcpp::List threshold_config, const Rcpp::List clust_config, const Rcpp::List parallel_config, const std::string output_type, const int verbose);
RcppExport SEXP _optimotu_seq_cluster_multi(SEXP seqSEXP, SEXP whichSEXP, SEXP dist_configSEXP, SEXP threshold_configSEXP, SEXP clust_configSEXP, SEXP parallel_configSEXP, SEXP output_typeSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::CharacterVector& >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< const Rcpp::ListOf<Rcpp::CharacterVector> >::type which(whichSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type dist_config(dist_configSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type threshold_config(threshold_configSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type clust_config(clust_configSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List >::type parallel_config(parallel_configSEXP);
    Rcpp::traits::input_parameter< const std::string >::type output_type(output_typeSEXP);
    Rcpp::traits::input_parameter< const int >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_cluster_multi(seq, which, dist_config, threshold_config, clust_config, parallel_config, output_type, verbose));
    return rcpp_result_gen;
END_RCPP
}
// seq_search_internal
Rcpp::RObject seq_search_internal(Rcpp::CharacterVector query, Rcpp::CharacterVector ref, Rcpp::List dist_config, Rcpp::List parallel_config, double threshold, int verbose, bool return_cigar, int span);
RcppExport SEXP _optimotu_seq_search_internal(SEXP querySEXP, SEXP refSEXP, SEXP dist_configSEXP, SEXP parallel_configSEXP, SEXP thresholdSEXP, SEXP verboseSEXP, SEXP return_cigarSEXP, SEXP spanSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type query(querySEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type ref(refSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type dist_config(dist_configSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type parallel_config(parallel_configSEXP);
    Rcpp::traits::input_parameter< double >::type threshold(thresholdSEXP);
    Rcpp::traits::input_parameter< int >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type return_cigar(return_cigarSEXP);
    Rcpp::traits::input_parameter< int >::type span(spanSEXP);
    rcpp_result_gen = Rcpp::wrap(seq_search_internal(query, ref, dist_config, parallel_config, threshold, verbose, return_cigar, span));
    return rcpp_result_gen;
END_RCPP
}
// summarize_by_rank
Rcpp::RObject summarize_by_rank(Rcpp::DataFrame data, Rcpp::CharacterVector ranks);
RcppExport SEXP _optimotu_summarize_by_rank(SEXP dataSEXP, SEXP ranksSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::DataFrame >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type ranks(ranksSEXP);
    rcpp_result_gen = Rcpp::wrap(summarize_by_rank(data, ranks));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_optimotu_confusion_matrix", (DL_FUNC) &_optimotu_confusion_matrix, 3},
    {"_optimotu_confusion_matrix2", (DL_FUNC) &_optimotu_confusion_matrix2, 3},
    {"_optimotu_distmx_cluster_single", (DL_FUNC) &_optimotu_distmx_cluster_single, 8},
    {"_optimotu_distmx_cluster_multi", (DL_FUNC) &_optimotu_distmx_cluster_multi, 9},
    {"_optimotu_seq_distmx_edlib", (DL_FUNC) &_optimotu_seq_distmx_edlib, 4},
    {"_optimotu_fastq_names", (DL_FUNC) &_optimotu_fastq_names, 1},
    {"_optimotu_intersect_length", (DL_FUNC) &_optimotu_intersect_length, 2},
    {"_optimotu_intersect_length_string", (DL_FUNC) &_optimotu_intersect_length_string, 2},
    {"_optimotu_inner_fmeasure", (DL_FUNC) &_optimotu_inner_fmeasure, 3},
    {"_optimotu_fmeasure_list", (DL_FUNC) &_optimotu_fmeasure_list, 3},
    {"_optimotu_fmeasure_matrix", (DL_FUNC) &_optimotu_fmeasure_matrix, 3},
    {"_optimotu_seq_distmx_kmer", (DL_FUNC) &_optimotu_seq_distmx_kmer, 10},
    {"_optimotu_mutual_information", (DL_FUNC) &_optimotu_mutual_information, 3},
    {"_optimotu_adjusted_mutual_information", (DL_FUNC) &_optimotu_adjusted_mutual_information, 3},
    {"_optimotu_pairwise_alignment", (DL_FUNC) &_optimotu_pairwise_alignment, 4},
    {"_optimotu_cigar_wfa2_global", (DL_FUNC) &_optimotu_cigar_wfa2_global, 8},
    {"_optimotu_cigar_wfa2_extend", (DL_FUNC) &_optimotu_cigar_wfa2_extend, 8},
    {"_optimotu_cigar_edlib_global", (DL_FUNC) &_optimotu_cigar_edlib_global, 2},
    {"_optimotu_cigar_edlib_extend", (DL_FUNC) &_optimotu_cigar_edlib_extend, 2},
    {"_optimotu_align_wfa2_global", (DL_FUNC) &_optimotu_align_wfa2_global, 8},
    {"_optimotu_align_wfa2_extend", (DL_FUNC) &_optimotu_align_wfa2_extend, 8},
    {"_optimotu_align_edlib_global", (DL_FUNC) &_optimotu_align_edlib_global, 2},
    {"_optimotu_align_edlib_extend", (DL_FUNC) &_optimotu_align_edlib_extend, 2},
    {"_optimotu_seq_distmx_wfa2", (DL_FUNC) &_optimotu_seq_distmx_wfa2, 11},
    {"_optimotu_seq_cluster_single", (DL_FUNC) &_optimotu_seq_cluster_single, 7},
    {"_optimotu_seq_cluster_multi", (DL_FUNC) &_optimotu_seq_cluster_multi, 8},
    {"_optimotu_seq_search_internal", (DL_FUNC) &_optimotu_seq_search_internal, 8},
    {"_optimotu_summarize_by_rank", (DL_FUNC) &_optimotu_summarize_by_rank, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_optimotu(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
