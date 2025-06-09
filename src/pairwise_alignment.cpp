#include "optimotu.h"
#include "pairwise_alignment.h"
#include "config.h"
#include <cstdint>
#include <sstream>
#include <string>
#include <optional>


// stream insertion operator for std::optional
// treats the stream as still valid even if nothing was found for the optional.
template<typename T>
std::istream& operator>>(std::istream& is, std::optional<T>& obj)
{
  if (T result; is >> result) {
    obj = result;
  }
  else {
    obj = {};
    is.clear();
  }
  return is;
}

// calculate distance from a cigar that may or may not be compressed.
double distance_from_cigar(const std::string & cigar) {
  if (cigar == "=") return 0.0;

  std::uint16_t match = 0, length = 0;
  std::optional<std::uint16_t> n;
  char op;
  std::istringstream ss(cigar);
  while (ss >> n >> op) {
    if (op == '=' || op == 'M') {
      match += n.value_or(1);
    }
    length += n.value_or(1);
  }
  return double(length - match) / double(length);
}

double distance_from_cigar_extend(const std::string & cigar) {
  if (cigar == "=") return 0.0;

  std::uint16_t match = 0, length = 0, end_gap = 0;
  std::optional<std::uint16_t> n;
  char op;
  std::istringstream ss(cigar);
  while (ss >> n >> op) {
    if (op == '=' || op == 'M') {
      match += n.value_or(1);
      end_gap = 0;
    } else if (op == 'I' || op == 'D') {
      end_gap += n.value_or(1);
    }
    length += n.value_or(1);
  }
  return double(length - match - end_gap) / double(length - end_gap);
}

template<enum AlignmentSpan span>
double distance_wfa2(const std::string &a, const std::string &b, wfa::WFAligner &aligner) {
  wfa::WFAligner::AlignmentStatus status;

  if constexpr (span == AlignmentSpan::GLOBAL) {
    status = aligner.alignEnd2End(a, b);
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    status = aligner.alignExtension(a, b);
  } else {
    static_assert(span != span, "Instantiation of unimplemented AlignmentSpan");
  }
  if (status != wfa::WFAligner::StatusAlgCompleted &&
      status != wfa::WFAligner::StatusAlgPartial) return 1.0;
  std::string cigar(aligner.getCIGAR(true));
  if constexpr (span == AlignmentSpan::EXTEND) {
    return distance_from_cigar_extend(cigar);
  } else {
    return distance_from_cigar(cigar);
  }
}

// explicit instantiations
template double distance_wfa2<AlignmentSpan::EXTEND>(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner
);

template double distance_wfa2<AlignmentSpan::GLOBAL>(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner
);

// calculate distance using edlib, with a pre-initialized config object
// this is
double distance_edlib(const std::string &a, const std::string &b, EdlibAlignConfig &aligner) {
  auto aln = edlibAlign(a.c_str(), a.size(), b.c_str(), b.size(), aligner);
  if (aln.status != EDLIB_STATUS_OK) return 1.0;
  if (aln.editDistance == -1) return 1.0;
  double d = (double)aln.editDistance / (double)aln.alignmentLength;
  edlibFreeAlignResult(aln);
  return d;
}

template<enum AlignmentSpan span>
double distance_edlib(const std::string &a, const std::string &b) {
  EdlibAlignMode mode;
  if constexpr (span == AlignmentSpan::GLOBAL) {
    mode = EDLIB_MODE_NW;
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    mode = EDLIB_MODE_SHW;
  } else {
    static_assert(span != span, "Instantiation of unimplemented AlignmentSpan");
  }
  auto aligner = edlibNewAlignConfig(-1, mode, EDLIB_TASK_PATH, NULL, 0);
  return distance_edlib(a, b, aligner);
}

//explicit instantiations
template double distance_edlib<AlignmentSpan::GLOBAL>(
    const std::string &a, const std::string &b);
template double distance_edlib<AlignmentSpan::EXTEND>(
    const std::string &a, const std::string &b);

template<enum AlignmentSpan span>
std::string cigar_wfa2(const std::string &a, const std::string &b,
                       wfa::WFAligner &aligner) {
  wfa::WFAligner::AlignmentStatus status;
  if constexpr (span == AlignmentSpan::GLOBAL) {
    status = aligner.alignEnd2End(a, b);
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    status = aligner.alignExtension(a, b);
  } else {
    static_assert(span != span, "Instantiation of unimplemented AlignmentSpan");
  }
  if (status != wfa::WFAligner::StatusAlgCompleted) return "";
  return aligner.getCIGAR(true);
}

std::string cigar_wfa2_global(const std::string &a, const std::string &b,
                              int match, int mismatch,
                              int gap_open, int gap_extend,
                              int gap_open2, int gap_extend2) {
  wfa::WFAlignerChoose aligner{match, mismatch, gap_open, gap_extend,
                               gap_open2, gap_extend2,
                               wfa::WFAligner::AlignmentScope::Alignment};
  return cigar_wfa2<AlignmentSpan::GLOBAL>(a, b, aligner);
}

std::string cigar_wfa2_extend(const std::string &a, const std::string &b,
                              int match, int mismatch,
                              int gap_open, int gap_extend,
                              int gap_open2, int gap_extend2) {
  wfa::WFAlignerChoose aligner{match, mismatch, gap_open, gap_extend,
                               gap_open2, gap_extend2,
                               wfa::WFAligner::AlignmentScope::Alignment};
  return cigar_wfa2<AlignmentSpan::EXTEND>(a, b, aligner);
}

template<enum AlignmentSpan span>
std::pair<double, std::string> distance_and_cigar_wfa2(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner
) {
  wfa::WFAligner::AlignmentStatus status;
  if constexpr (span == AlignmentSpan::GLOBAL) {
    status = aligner.alignEnd2End(a, b);
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    status = aligner.alignExtension(a, b);
  } else {
    static_assert(span != span, "Instantiation of unimplemented AlignmentSpan");
  }

  if (status != wfa::WFAligner::StatusAlgCompleted) return {1.0, ""};
  std::string cigar(aligner.getCIGAR(true));
  if constexpr (span == AlignmentSpan::EXTEND) {
    return {distance_from_cigar_extend(cigar), cigar};
  } else {
    return {distance_from_cigar(cigar), cigar};
  }
}

// explicit instantiations
template
std::pair<double, std::string> distance_and_cigar_wfa2<AlignmentSpan::EXTEND>(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner
);

template
std::pair<double, std::string> distance_and_cigar_wfa2<AlignmentSpan::GLOBAL>(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner
);

template<enum AlignmentSpan span>
std::string cigar_edlib(const std::string &a, const std::string &b) {

  EdlibAlignMode mode;
  if constexpr (span == AlignmentSpan::GLOBAL) {
    mode = EDLIB_MODE_NW;
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    mode = EDLIB_MODE_SHW;
  } else {
    static_assert(span != span, "Instantiation of unimplemented AlignmentSpan");
  }
  auto config = edlibNewAlignConfig(-1, mode, EDLIB_TASK_PATH, NULL, 0);

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

std::string cigar_edlib_global(const std::string &a, const std::string &b) {
  return cigar_edlib<AlignmentSpan::GLOBAL>(a, b);
}

std::string cigar_edlib_extend(const std::string &a, const std::string &b) {
  return cigar_edlib<AlignmentSpan::EXTEND>(a, b);
}

template<enum AlignmentSpan span>
std::pair<int, double> score_and_distance_wfa2(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner
) {
  wfa::WFAligner::AlignmentStatus status;
  // if (a.size() >= b.size()) {
  if constexpr (span == AlignmentSpan::GLOBAL) {
    status = aligner.alignEnd2End(a, b);
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    status = aligner.alignExtension(a, b);
  } else {
    static_assert(span != span, "Instantiation of unimplemented AlignmentSpan");
  }
  // } else {
  //   status = aligner.alignEnd2End(b, a);
  // }
  if (status != wfa::WFAligner::StatusAlgCompleted) return {
    std::max(a.size(), b.size()), 1.0};
  auto cigar = aligner.getCIGAR(true);
  if constexpr (span == AlignmentSpan::EXTEND) {
    return {aligner.getAlignmentScore(), distance_from_cigar_extend(cigar)};
  } else {
    return {aligner.getAlignmentScore(), distance_from_cigar(cigar)};
  }
}

// explicit instantiations
template std::pair<int, double> score_and_distance_wfa2<AlignmentSpan::GLOBAL>(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner
);

template std::pair<int, double> score_and_distance_wfa2<AlignmentSpan::EXTEND>(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner
);


double align_wfa2_global(const std::string a, const std::string b,
             int match, int mismatch,
             int gap_open, int gap_extend,
             int gap_open2, int gap_extend2) {
  wfa::WFAlignerChoose aligner{match, mismatch, gap_open, gap_extend,
                               gap_open2, gap_extend2,
                               wfa::WFAligner::AlignmentScope::Alignment};
  return distance_wfa2<AlignmentSpan::GLOBAL>(a, b, aligner);
}

double align_wfa2_extend(const std::string a, const std::string b,
             int match, int mismatch,
             int gap_open, int gap_extend,
             int gap_open2, int gap_extend2) {
  wfa::WFAlignerChoose aligner{match, mismatch, gap_open, gap_extend,
                               gap_open2, gap_extend2,
                               wfa::WFAligner::AlignmentScope::Alignment};
  return distance_wfa2<AlignmentSpan::EXTEND>(a, b, aligner);
}

double align_edlib_global(const std::string a, const std::string b) {
  auto config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);
  return distance_edlib(a, b, config);
}

double align_edlib_extend(const std::string a, const std::string b) {
  auto config = edlibNewAlignConfig(-1, EDLIB_MODE_SHW, EDLIB_TASK_PATH, NULL, 0);
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
  out[0].reserve(std::max({cigar.size(), a.size(), b.size()}));
  out[1].reserve(out[0].capacity());
  size_t i = 0, j = 0;
  size_t k = 0;
  while (k < cigar.size()) {
    // if there are digits, parse them
    size_t n = 0;
    while (k < cigar.size() && std::isdigit(cigar[k])) {
      n = 10 * n + (cigar[k] - '0');
      k++;
    }
    // if no digits, assume 1
    if (n == 0) {
      n = 1;
    }
    if (k >= cigar.size()) {
      OPTIMOTU_STOP("CIGAR string ends with a number");
    }
    switch (cigar[k]) {
    case '=':
    case 'M':
    case 'X':
      for (size_t m = 0; m < n; m++) {
        out[0].push_back(a[i]);
        out[1].push_back(b[j]);
        i++;
        j++;
      }
      break;
    case 'D':
      for (size_t m = 0; m < n; m++) {
        out[0].push_back(a[i]);
        out[1].push_back('-');
        i++;
      }
      break;
    case 'I':
      for (size_t m = 0; m < n; m++) {
        out[0].push_back('-');
        out[1].push_back(b[j]);
        j++;
      }
      break;
    default:
      OPTIMOTU_STOP("Unknown CIGAR character: " + std::string(1, cigar[k]));
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
  Rcpp::List dist_config,
  int span = 0
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
    if (span == 0) {
      return align_from_compressed_cigar(a, b, cigar_wfa2_global(a, b, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2));
    } else if (span == 1) {
      return align_from_compressed_cigar(a, b, cigar_wfa2_extend(a, b, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2));
    } else {
      OPTIMOTU_STOP(
        "Unknown alignment span: '" + std::to_string(span) + "'"
      );
    }
  } else if (dist_method == "edlib") {
    if (span == 0) {
      return align_from_compressed_cigar(a, b, cigar_edlib_global(a, b));
    } else if (span == 1) {
      return align_from_compressed_cigar(a, b, cigar_edlib_extend(a, b));
    } else {
      OPTIMOTU_STOP(
        "Unknown alignment span: '" + std::to_string(span) + "'"
      );
    }
  } else {
    OPTIMOTU_STOP(
      "Unknown distance method: '" + dist_method + "'"
    );
  }
}
