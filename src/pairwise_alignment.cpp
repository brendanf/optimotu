#include "pairwise_alignment.h"
#include <cstdint>

double distance_wfa2(const std::string &a, const std::string &b, wfa::WFAligner &aligner) {
  auto status = aligner.alignEnd2End(a, b);
  if (status != wfa::WFAligner::StatusAlgCompleted) return 1.0;
  auto cigar = aligner.getCIGAR(true);
  std::uint16_t match = 0, length = 0;
  for (char c : cigar) {
    if (c == 'M') {
      match++;
    }
    length++;
  }
  return 1.0 - double(match) / double(length);
}

double distance_wfa2(const std::string&a, const std::string&b, wfa::WFAlignerEdit &aligner) {
  auto status = aligner.alignEnd2End(a, b);
  if (status != wfa::WFAligner::StatusAlgCompleted) return 1.0;
  auto cigar = aligner.getCIGAR(true);
  return double(aligner.getAlignmentScore()) / double(cigar.size());
}

double distance_edlib(const std::string &a, const std::string &b, EdlibAlignConfig &aligner) {
  auto aln = edlibAlign(a.c_str(), a.size(), b.c_str(), b.size(), aligner);
  if (aln.status != EDLIB_STATUS_OK) return 1.0;
  if (aln.editDistance == -1) return 1.0;
  double d = (double)aln.editDistance / (double)aln.alignmentLength;
  edlibFreeAlignResult(aln);
  return d;
}

std::string cigar_wfa2(const std::string &a, const std::string &b,
                       int match, int mismatch,
                       int open1, int extend1,
                       int open2, int extend2) {
  wfa::WFAlignerChoose aligner{match, mismatch, open1, extend1, open2, extend2,
                               wfa::WFAligner::AlignmentScope::Alignment};
  wfa::WFAligner::AlignmentStatus status;
  // if (a.size() >= b.size()) {
  status = aligner.alignEnd2End(a, b);
  // } else {
  //   status = aligner.alignEnd2End(b, a);
  // }
  switch (status) {
  case wfa::WFAligner::StatusOOM:
#ifdef OPTIMOTU_R
    Rcpp::stop("Aligner out of memory.");
#else
    std::cout << "Aligner out of memory." << std::endl;
    exit(1);
#endif
    break;
  case wfa::WFAligner::StatusAlgPartial:
  case wfa::WFAligner::StatusMaxStepsReached:
#ifdef OPTIMOTU_R
    Rcpp::stop("Alignment not feasible with given constraints.");
#else
    std::cout << "Alignment not feasible with given constraints." << std::endl;
    exit(1);
#endif
    break;
  case wfa::WFAligner::StatusAlgCompleted:
    break;
  }
  return aligner.getCIGAR(true);
}

std::string cigar_edlib(const std::string &a, const std::string &b) {
  auto config = edlibNewAlignConfig(-1, EdlibAlignMode::EDLIB_MODE_NW,
                                    EdlibAlignTask::EDLIB_TASK_PATH, NULL, 0);
  std::vector<char> key;
  EdlibAlignResult alignment;
  // if (a.size() >= b.size()){
  alignment = edlibAlign(a.c_str(), a.size(), b.c_str(), b.size(),
                         config);
  key = {'M', 'D', 'I', 'X'};
  // } else {
  //   alignment = edlibAlign(b.c_str(), b.size(), a.c_str(), a.size(),
  //                          config);
  //   key = {'M', 'I', 'D', 'X'};
  //
  // }
  // Rcpp::Rcout << "a_len=" << a.size() << std::endl
  //             << "b_len=" << b.size() << std::endl
  //             << "aln_len=" << alignment.alignmentLength << std::endl
  //             << "start_pos=" << alignment.startLocations[0] << std::endl
  //             << "end_pos=" << alignment.endLocations[0] << std::endl;

  std::string out{};
  out.reserve(alignment.alignmentLength);
  for (int i = 0; i < alignment.alignmentLength; ++i) {
    out.push_back(key[alignment.alignment[i]]);
  }
  edlibFreeAlignResult(alignment);
  return out;
}

std::pair<int, double> score_and_distance_wfa2(const std::string &a, const std::string &b, wfa::WFAligner &aligner) {
  wfa::WFAligner::AlignmentStatus status;
  // if (a.size() >= b.size()) {
  status = aligner.alignEnd2End(a, b);
  // } else {
  //   status = aligner.alignEnd2End(b, a);
  // }
  if (status != wfa::WFAligner::StatusAlgCompleted) return {std::max(a.size(), b.size()), 1.0};
  auto cigar = aligner.getCIGAR(true);
  std::uint16_t match = 0, length = 0;
  for (char c : cigar) {
    if (c == 'M') {
      match++;
    }
    length++;
  }
  return {length - match, double(length - match) / double(length)};
}

std::pair<int, double> score_and_distance_wfa2(const std::string &a, const std::string &b, wfa::WFAlignerEdit &aligner) {
  wfa::WFAligner::AlignmentStatus status;
  // if (a.size() >= b.size()) {
  status = aligner.alignEnd2End(a, b);
  // } else {
  //   status = aligner.alignEnd2End(b, a);
  // }
  if (status != wfa::WFAligner::StatusAlgCompleted) return {std::max(a.size(), b.size()), 1.0};
  auto cigar = aligner.getCIGAR(true);
  return {aligner.getAlignmentScore(),
          double(aligner.getAlignmentScore()) / double(cigar.size())};
}

double align(const std::string a, const std::string b,
             int match, int mismatch,
             int gap, int extend,
             int gap2, int extend2) {
  if (gap2 != 0 || extend2 != 0) {
    wfa::WFAlignerGapAffine2Pieces aligner(
        match, mismatch,
        gap, extend,
        gap2, extend2,
        wfa::WFAligner::Alignment);
    return distance_wfa2(a, b, aligner);
  } else if (extend != 0) {
    wfa::WFAlignerGapAffine aligner(
        match, mismatch,
        gap, extend,
        wfa::WFAligner::Alignment);
    return distance_wfa2(a, b, aligner);
  } else if (match == 0 && gap == mismatch) {
    wfa::WFAlignerEdit aligner(wfa::WFAligner::Alignment);
    return distance_wfa2(a, b, aligner);
  } else if (mismatch == 0 && match == 0) {
    wfa::WFAlignerIndel aligner(wfa::WFAligner::Alignment);
    return distance_wfa2(a, b, aligner);
  } else {
    wfa::WFAlignerGapLinear aligner(
        match, mismatch, gap,
        wfa::WFAligner::Alignment);
    return distance_wfa2(a, b, aligner);
  }
}
