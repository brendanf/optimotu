#include "optimotu.h"
#include "pairwise_alignment.h"
#include "config.h"
#include <cstdint>
#include <sstream>
#include <string>

double distance_wfa2(const std::string &a, const std::string &b, wfa::WFAligner &aligner) {
  auto status = aligner.alignEnd2End(a, b);
  if (status != wfa::WFAligner::StatusAlgCompleted) return 1.0;
  std::string cigar(aligner.getCIGAR(true));
  std::uint16_t match = 0, length = 0;
  std::uint16_t l;
  char c;
  std::istringstream ss(cigar);
  while (ss >> l >> c) {
    if (c == '=') {
      match += l;
    }
    length += l;
  }
  return 1.0 - double(match) / double(length);
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
                       int gap_open, int gap_extend,
                       int gap_open2, int gap_extend2) {
  wfa::WFAlignerChoose aligner{-match, mismatch, gap_open, gap_extend,
                               gap_open2, gap_extend2,
                               wfa::WFAligner::AlignmentScope::Alignment};
  wfa::WFAligner::AlignmentStatus status;
  // if (a.size() >= b.size()) {
  status = aligner.alignEnd2End(a, b);
  // } else {
  //   status = aligner.alignEnd2End(b, a);
  // }
  switch (status) {
  case wfa::WFAligner::StatusOOM:
    OPTIMOTU_STOP("Aligner out of memory.");
    break;
  case wfa::WFAligner::StatusAlgPartial:
  case wfa::WFAligner::StatusMaxStepsReached:
    OPTIMOTU_STOP("Alignment not feasible with given constraints.");
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

double align_wfa2(const std::string a, const std::string b,
             int match, int mismatch,
             int gap_open, int gap_extend,
             int gap_open2, int gap_extend2) {
  wfa::WFAlignerChoose aligner{match, mismatch, gap_open, gap_extend,
                               gap_open2, gap_extend2,
                               wfa::WFAligner::AlignmentScope::Alignment};
  return distance_wfa2(a, b, aligner);
}

double align_edlib(const std::string a, const std::string b) {
  auto config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);
  return distance_edlib(a, b, config);
}

std::vector<std::string> align_from_cigar(
    const std::string & a,
    const std::string & b,
    const std::string & cigar
) {
  std::vector<std::string> out(2);
  out[0].reserve(cigar.size());
  out[1].reserve(cigar.size());
  size_t i = 0, j = 0;
  for (char c : cigar) {
    switch (c) {
    case 'M':
      out[0].push_back(a[i]);
      out[1].push_back(b[j]);
      i++;
      j++;
      break;
    case 'D':
      out[0].push_back(a[i]);
      out[1].push_back('-');
      i++;
      break;
    case 'I':
      out[0].push_back('-');
      out[1].push_back(b[j]);
      j++;
      break;
    case 'X':
      out[0].push_back(a[i]);
      out[1].push_back(b[j]);
      i++;
      j++;
      break;
    }
  }
  return out;
}

std::vector<std::string> align_from_compressed_cigar(
    const std::string & a,
    const std::string & b,
    const std::string & cigar
) {
  std::vector<std::string> out(2);
  out[0].reserve(cigar.size());
  out[1].reserve(cigar.size());
  size_t i = 0, j = 0;
  size_t k = 0;
  while (k < cigar.size()) {
    size_t l = 0;
    while (k + l < cigar.size() && cigar[k + l] >= '0' && cigar[k + l] <= '9') {
      l++;
    }
    int n = std::stoi(cigar.substr(k, l));
    k += l;
    switch (cigar[k]) {
    case '=':
      for (int m = 0; m < n; m++) {
        out[0].push_back(a[i]);
        out[1].push_back(b[j]);
        i++;
        j++;
      }
      break;
    case 'D':
      for (int m = 0; m < n; m++) {
        out[0].push_back(a[i]);
        out[1].push_back('-');
        i++;
      }
      break;
    case 'I':
      for (int m = 0; m < n; m++) {
        out[0].push_back('-');
        out[1].push_back(b[j]);
        j++;
      }
      break;
    case 'X':
      for (int m = 0; m < n; m++) {
        out[0].push_back(a[i]);
        out[1].push_back(b[j]);
        i++;
        j++;
      }
      break;
    }
    k++;
  }
  return out;

}

// [[Rcpp::export]]
std::vector<std::string> pairwise_alignment(
  std::string a,
  std::string b,
  Rcpp::List dist_config
) {
  if (!dist_config.inherits("optimotu_dist_config")) {
    OPTIMOTU_STOP(
      "'dist_config' must be of class 'optimotu_dist_config'"
    );
  }
  std::string dist_method = element_as_string(dist_config, "method", "dist_config");
  if (dist_method == "wfa2") {
    int match = element_as_int(dist_config, "match", "dist_config");
    int mismatch = element_as_int(dist_config, "mismatch", "dist_config");
    int gap_open = element_as_int(dist_config, "gap_open", "dist_config");
    int gap_extend = element_as_int(dist_config, "gap_extend", "dist_config");
    int gap_open2 = element_as_int(dist_config, "gap_open2", "dist_config");
    int gap_extend2 = element_as_int(dist_config, "gap_extend2", "dist_config");
    return align_from_compressed_cigar(a, b, cigar_wfa2(a, b, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2));
  } else if (dist_method == "edlib") {
    return align_from_cigar(a, b, cigar_edlib(a, b));
  } else {
    OPTIMOTU_STOP(
      "Unknown distance method: '" + dist_method + "'"
    );
  }
}
