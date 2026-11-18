#include "optimotu.h"
#include "pairwise_alignment.h"
#include "config.h"
#include "score_dissimilarity.h"
#include "wfa_fatal.h"
#include <cstdint>
#include <sstream>
#include <string>
#include <optional>

namespace {

Dissimilarity parse_dissimilarity(const std::string &s) {
  if (s == "score") return Dissimilarity::SCORE;
  if (s == "identity") return Dissimilarity::IDENTITY;
  OPTIMOTU_STOP("dissimilarity must be \"identity\" or \"score\"");
}

bool wfa_hard_fail(
    wfa::WFAligner &aligner,
    wfa::WFAligner::AlignmentStatus status)
{
  if (aligner.fatal() || status == wfa::WFAligner::StatusOOM)
  {
    optimotu_wfa_note_fatal();
    return true;
  }
  return false;
}

}  // namespace


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

std::tuple<double, int, int, int, int, int> dist_gapstats_from_cigar(const std::string & cigar) {
  // difficult case because the CIGAR string does not tell the length
  if (cigar == "=") return {0.0, 0, 0, 0, 0, 0};
  std::uint16_t match = 0;
  int length = 0, n_insert = 0, n_delete = 0, max_insert = 0, max_delete = 0, ins_run = 0, del_run = 0;
  std::optional<std::uint16_t> n;
  char op;
  std::istringstream ss(cigar);
  while (ss >> n >> op) {
    if (op == '=' || op == 'M') {
      match += n.value_or(1);
      ins_run = 0;
      del_run = 0;
    } else if (op == 'I') {
      n_insert += n.value_or(1);
      ins_run += n.value_or(1);
      if (ins_run > max_insert) max_insert = ins_run;
    } else if (op == 'D') {
      n_delete += n.value_or(1);
      del_run += n.value_or(1);
      if (del_run > max_delete) max_delete = del_run;
    } else {
      ins_run = 0;
      del_run = 0;
    }
    length += n.value_or(1);
  }
  return {double(length - match) / double(length), length, n_insert, n_delete, max_insert, max_delete};
}

std::tuple<double, int, int, int, int, int> dist_gapstats_from_cigar_extend(const std::string & cigar) {
  if (cigar == "=") return {0.0, 0, 0, 0, 0, 0};
  std::uint16_t match = 0, end_ins = 0, end_del = 0;
  int length = 0, n_insert = 0, n_delete = 0, max_insert = 0, max_delete = 0, ins_run = 0, del_run = 0;
  std::optional<std::uint16_t> n;
  char op;
  std::istringstream ss(cigar);
  while (ss >> n >> op) {
    if (op == '=' || op == 'M') {
      match += n.value_or(1);
      if (ins_run > max_insert) max_insert = ins_run;
      if (del_run > max_delete) max_delete = del_run;
      ins_run = 0;
      del_run = 0;
      end_ins = 0;
      end_del = 0;
    } else if (op == 'I') {
      end_ins += n.value_or(1);
      n_insert += n.value_or(1);
      ins_run += n.value_or(1);
    } else if (op == 'D') {
      end_del += n.value_or(1);
      n_delete += n.value_or(1);
      del_run += n.value_or(1);
    } else {
      if (ins_run > max_insert) max_insert = ins_run;
      if (del_run > max_delete) max_delete = del_run;
      ins_run = 0;
      del_run = 0;
      end_ins = 0;
      end_del = 0;
    }
    length += n.value_or(1);
  }
  return {double(length - match - end_ins - end_del) / double(length - end_ins - end_del), length, n_insert - end_ins, n_delete - end_del, max_insert, max_delete};
}


template<enum AlignmentSpan span>
double distance_wfa2(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner,
    Dissimilarity dissimilarity,
    int match,
    int mismatch
) {
  wfa::WFAligner::AlignmentStatus status;

  if constexpr (span == AlignmentSpan::GLOBAL) {
    status = aligner.alignEnd2End(a, b);
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    status = aligner.alignExtension(a, b);
  } else {
    static_assert(always_false_v<span>, "Instantiation of unimplemented AlignmentSpan");
  }
  if (wfa_hard_fail(aligner, status))
    return 1.0;
  if (status != wfa::WFAligner::StatusAlgCompleted &&
      status != wfa::WFAligner::StatusAlgPartial) return 1.0;
  if (dissimilarity == Dissimilarity::SCORE) {
    int L = score_ref_length(a.size(), b.size(), span);
    int cost = wfa_reported_to_cost(aligner.getAlignmentScore(), match);
    return cost_score_distance(cost, L, match, mismatch);
  }
  std::string cigar(aligner.getCIGAR(true));
  if (aligner.fatal())
  {
    optimotu_wfa_note_fatal();
    return 1.0;
  }
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
    wfa::WFAligner &aligner,
    Dissimilarity dissimilarity,
    int match,
    int mismatch
);

template double distance_wfa2<AlignmentSpan::GLOBAL>(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner,
    Dissimilarity dissimilarity,
    int match,
    int mismatch
);

template<enum AlignmentSpan span>
double distance_wfa2(
    const SequenceView &a,
    const SequenceView &b,
    wfa::WFAligner &aligner,
    Dissimilarity dissimilarity,
    int match,
    int mismatch
) {
  wfa::WFAligner::AlignmentStatus status;

  if constexpr (span == AlignmentSpan::GLOBAL) {
    status = aligner.alignEnd2End(
      a.data(), static_cast<int>(a.size()),
      b.data(), static_cast<int>(b.size())
    );
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    status = aligner.alignExtension(
      a.data(), static_cast<int>(a.size()),
      b.data(), static_cast<int>(b.size())
    );
  } else {
    static_assert(always_false_v<span>, "Instantiation of unimplemented AlignmentSpan");
  }
  if (wfa_hard_fail(aligner, status))
    return 1.0;
  if (status != wfa::WFAligner::StatusAlgCompleted &&
      status != wfa::WFAligner::StatusAlgPartial) return 1.0;
  if (dissimilarity == Dissimilarity::SCORE) {
    int L = score_ref_length(a.size(), b.size(), span);
    int cost = wfa_reported_to_cost(aligner.getAlignmentScore(), match);
    return cost_score_distance(cost, L, match, mismatch);
  }
  std::string cigar(aligner.getCIGAR(true));
  if (aligner.fatal())
  {
    optimotu_wfa_note_fatal();
    return 1.0;
  }
  if constexpr (span == AlignmentSpan::EXTEND) {
    return distance_from_cigar_extend(cigar);
  } else {
    return distance_from_cigar(cigar);
  }
}

template double distance_wfa2<AlignmentSpan::EXTEND>(
    const SequenceView &a,
    const SequenceView &b,
    wfa::WFAligner &aligner,
    Dissimilarity dissimilarity,
    int match,
    int mismatch
);

template double distance_wfa2<AlignmentSpan::GLOBAL>(
    const SequenceView &a,
    const SequenceView &b,
    wfa::WFAligner &aligner,
    Dissimilarity dissimilarity,
    int match,
    int mismatch
);

static double distance_from_edlib_result(
    const EdlibAlignResult &aln,
    int l1,
    int l2,
    EdlibAlignMode mode,
    Dissimilarity dissimilarity
) {
  if (aln.status != EDLIB_STATUS_OK || aln.editDistance < 0) return 1.0;
  if (dissimilarity == Dissimilarity::SCORE) {
    AlignmentSpan span = (mode == EDLIB_MODE_SHW)
      ? AlignmentSpan::EXTEND
      : AlignmentSpan::GLOBAL;
    int L = score_ref_length(l1, l2, span);
    return cost_score_distance(aln.editDistance, L, 0, 1);
  }
  if (aln.alignmentLength <= 0) return 1.0;
  return (double)aln.editDistance / (double)aln.alignmentLength;
}

double distance_edlib(
    const std::string &a,
    const std::string &b,
    EdlibAlignConfig &aligner,
    Dissimilarity dissimilarity
) {
  auto aln = edlibAlign(a.c_str(), a.size(), b.c_str(), b.size(), aligner);
  double d = distance_from_edlib_result(
      aln, static_cast<int>(a.size()), static_cast<int>(b.size()),
      aligner.mode, dissimilarity
  );
  edlibFreeAlignResult(aln);
  return d;
}

double distance_edlib(
    const SequenceView &a,
    const SequenceView &b,
    EdlibAlignConfig &aligner,
    Dissimilarity dissimilarity
) {
  auto aln = edlibAlign(
    a.data(), static_cast<int>(a.size()),
    b.data(), static_cast<int>(b.size()),
    aligner
  );
  double d = distance_from_edlib_result(
      aln, static_cast<int>(a.size()), static_cast<int>(b.size()),
      aligner.mode, dissimilarity
  );
  edlibFreeAlignResult(aln);
  return d;
}

template<enum AlignmentSpan span>
double distance_edlib(
    const std::string &a,
    const std::string &b,
    Dissimilarity dissimilarity
) {
  EdlibAlignMode mode;
  if constexpr (span == AlignmentSpan::GLOBAL) {
    mode = EDLIB_MODE_NW;
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    mode = EDLIB_MODE_SHW;
  } else {
    static_assert(always_false_v<span>, "Instantiation of unimplemented AlignmentSpan");
  }
  EdlibAlignTask task = (dissimilarity == Dissimilarity::SCORE)
    ? EDLIB_TASK_DISTANCE
    : EDLIB_TASK_PATH;
  auto aligner = edlibNewAlignConfig(-1, mode, task, NULL, 0);
  return distance_edlib(a, b, aligner, dissimilarity);
}

//explicit instantiations
template double distance_edlib<AlignmentSpan::GLOBAL>(
    const std::string &a, const std::string &b, Dissimilarity dissimilarity);
template double distance_edlib<AlignmentSpan::EXTEND>(
    const std::string &a, const std::string &b, Dissimilarity dissimilarity);

template<enum AlignmentSpan span>
std::string cigar_wfa2(const std::string &a, const std::string &b,
                       wfa::WFAligner &aligner) {
  wfa::WFAligner::AlignmentStatus status;
  if constexpr (span == AlignmentSpan::GLOBAL) {
    status = aligner.alignEnd2End(a, b);
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    status = aligner.alignExtension(a, b);
  } else {
    static_assert(always_false_v<span>, "Instantiation of unimplemented AlignmentSpan");
  }
  if (wfa_hard_fail(aligner, status))
    return "";
  if (status != wfa::WFAligner::StatusAlgCompleted) return "";
  std::string cigar = aligner.getCIGAR(true);
  if (aligner.fatal())
  {
    optimotu_wfa_note_fatal();
    return "";
  }
  return cigar;
}

std::string cigar_wfa2_global(const std::string &a, const std::string &b,
                              int match, int mismatch,
                              int gap_open, int gap_extend,
                              int gap_open2, int gap_extend2) {
  wfa::WFAlignerChoose aligner{match, mismatch, gap_open, gap_extend,
                               gap_open2, gap_extend2,
                               wfa::WFAligner::AlignmentScope::Alignment};
  std::string cigar = cigar_wfa2<AlignmentSpan::GLOBAL>(a, b, aligner);
  optimotu_wfa_stop_if_fatal();
  return cigar;
}

std::string cigar_wfa2_extend(const std::string &a, const std::string &b,
                              int match, int mismatch,
                              int gap_open, int gap_extend,
                              int gap_open2, int gap_extend2) {
  wfa::WFAlignerChoose aligner{match, mismatch, gap_open, gap_extend,
                               gap_open2, gap_extend2,
                               wfa::WFAligner::AlignmentScope::Alignment};
  std::string cigar = cigar_wfa2<AlignmentSpan::EXTEND>(a, b, aligner);
  optimotu_wfa_stop_if_fatal();
  return cigar;
}

// explicit instantiations
// shouldn't be needed, but without them I get linker errors
template std::string cigar_wfa2<AlignmentSpan::GLOBAL>(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner
);
template std::string cigar_wfa2<AlignmentSpan::EXTEND>(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner
);

template<enum AlignmentSpan span>
std::pair<double, std::string> distance_and_cigar_wfa2(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner,
    Dissimilarity dissimilarity,
    int match,
    int mismatch
) {
  wfa::WFAligner::AlignmentStatus status;
  if constexpr (span == AlignmentSpan::GLOBAL) {
    status = aligner.alignEnd2End(a, b);
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    status = aligner.alignExtension(a, b);
  } else {
    static_assert(always_false_v<span>, "Instantiation of unimplemented AlignmentSpan");
  }

  if (wfa_hard_fail(aligner, status))
    return {1.0, ""};
  if (status != wfa::WFAligner::StatusAlgCompleted) return {1.0, ""};
  std::string cigar(aligner.getCIGAR(true));
  if (aligner.fatal())
  {
    optimotu_wfa_note_fatal();
    return {1.0, ""};
  }
  double d;
  if (dissimilarity == Dissimilarity::SCORE) {
    int L = score_ref_length(a.size(), b.size(), span);
    int cost = wfa_reported_to_cost(aligner.getAlignmentScore(), match);
    d = cost_score_distance(cost, L, match, mismatch);
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    d = distance_from_cigar_extend(cigar);
  } else {
    d = distance_from_cigar(cigar);
  }
  return {d, cigar};
}

// explicit instantiations
template
std::pair<double, std::string> distance_and_cigar_wfa2<AlignmentSpan::EXTEND>(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner,
    Dissimilarity dissimilarity,
    int match,
    int mismatch
);

template
std::pair<double, std::string> distance_and_cigar_wfa2<AlignmentSpan::GLOBAL>(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner,
    Dissimilarity dissimilarity,
    int match,
    int mismatch
);

template<enum AlignmentSpan span>
std::string cigar_edlib(const std::string &a, const std::string &b) {

  EdlibAlignMode mode;
  if constexpr (span == AlignmentSpan::GLOBAL) {
    mode = EDLIB_MODE_NW;
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    mode = EDLIB_MODE_SHW;
  } else {
    static_assert(always_false_v<span>, "Instantiation of unimplemented AlignmentSpan");
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
std::pair<double, std::string> distance_and_cigar_edlib(
    const std::string &a,
    const std::string &b,
    EdlibAlignConfig &aligner,
    Dissimilarity dissimilarity
) {
  aligner.task = EDLIB_TASK_PATH;
  auto alignment = edlibAlign(a.c_str(), a.size(), b.c_str(), b.size(), aligner);
  if (alignment.status != EDLIB_STATUS_OK) return {1.0, ""};
  if (alignment.editDistance == -1) return {1.0, ""};
  double d = distance_from_edlib_result(
      alignment, static_cast<int>(a.size()), static_cast<int>(b.size()),
      aligner.mode, dissimilarity
  );
  std::string cigar;
  cigar.reserve(alignment.alignmentLength);
  std::vector<char> key = {'M', 'D', 'I', 'X'};
  for (int k = 0; k < alignment.alignmentLength; ++k) {
    cigar.push_back(key[alignment.alignment[k]]);
  }
  edlibFreeAlignResult(alignment);
  return {d, cigar};
}

// explicit instantiations
template std::pair<double, std::string> distance_and_cigar_edlib<AlignmentSpan::GLOBAL>(
    const std::string &a,
    const std::string &b,
    EdlibAlignConfig &aligner,
    Dissimilarity dissimilarity
);
template std::pair<double, std::string> distance_and_cigar_edlib<AlignmentSpan::EXTEND>(
    const std::string &a,
    const std::string &b,
    EdlibAlignConfig &aligner,
    Dissimilarity dissimilarity
);

template<enum AlignmentSpan span>
std::pair<int, double> score_and_distance_wfa2(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner,
    Dissimilarity dissimilarity,
    int match,
    int mismatch
) {
  wfa::WFAligner::AlignmentStatus status;
  if constexpr (span == AlignmentSpan::GLOBAL) {
    status = aligner.alignEnd2End(a, b);
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    status = aligner.alignExtension(a, b);
  } else {
    static_assert(always_false_v<span>, "Instantiation of unimplemented AlignmentSpan");
  }
  if (wfa_hard_fail(aligner, status))
    return {
        static_cast<int>(std::max(a.size(), b.size())), 1.0};
  if (status != wfa::WFAligner::StatusAlgCompleted) return {
    static_cast<int>(std::max(a.size(), b.size())), 1.0};
  if (aligner.fatal())
  {
    optimotu_wfa_note_fatal();
    return {static_cast<int>(std::max(a.size(), b.size())), 1.0};
  }
  int score = aligner.getAlignmentScore();
  if (dissimilarity == Dissimilarity::SCORE) {
    int L = score_ref_length(a.size(), b.size(), span);
    int cost = wfa_reported_to_cost(score, match);
    return {cost, cost_score_distance(cost, L, match, mismatch)};
  }
  auto cigar = aligner.getCIGAR(true);
  if constexpr (span == AlignmentSpan::EXTEND) {
    return {score, distance_from_cigar_extend(cigar)};
  } else {
    return {score, distance_from_cigar(cigar)};
  }
}

// explicit instantiations
template std::pair<int, double> score_and_distance_wfa2<AlignmentSpan::GLOBAL>(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner,
    Dissimilarity dissimilarity,
    int match,
    int mismatch
);

template std::pair<int, double> score_and_distance_wfa2<AlignmentSpan::EXTEND>(
    const std::string &a,
    const std::string &b,
    wfa::WFAligner &aligner,
    Dissimilarity dissimilarity,
    int match,
    int mismatch
);


double align_wfa2_global(const std::string a, const std::string b,
             int match, int mismatch,
             int gap_open, int gap_extend,
             int gap_open2, int gap_extend2,
             std::string dissimilarity) {
  Dissimilarity mode = parse_dissimilarity(dissimilarity);
  auto scope = mode == Dissimilarity::SCORE
    ? wfa::WFAligner::AlignmentScope::Score
    : wfa::WFAligner::AlignmentScope::Alignment;
  wfa::WFAlignerChoose aligner{match, mismatch, gap_open, gap_extend,
                               gap_open2, gap_extend2, scope};
  double d = distance_wfa2<AlignmentSpan::GLOBAL>(
      a, b, aligner, mode, match, mismatch);
  optimotu_wfa_stop_if_fatal();
  return d;
}

double align_wfa2_extend(const std::string a, const std::string b,
             int match, int mismatch,
             int gap_open, int gap_extend,
             int gap_open2, int gap_extend2,
             std::string dissimilarity) {
  Dissimilarity mode = parse_dissimilarity(dissimilarity);
  auto scope = mode == Dissimilarity::SCORE
    ? wfa::WFAligner::AlignmentScope::Score
    : wfa::WFAligner::AlignmentScope::Alignment;
  wfa::WFAlignerChoose aligner{match, mismatch, gap_open, gap_extend,
                               gap_open2, gap_extend2, scope};
  double d = distance_wfa2<AlignmentSpan::EXTEND>(
      a, b, aligner, mode, match, mismatch);
  optimotu_wfa_stop_if_fatal();
  return d;
}

double align_edlib_global(
    const std::string a,
    const std::string b,
    std::string dissimilarity
) {
  Dissimilarity mode = parse_dissimilarity(dissimilarity);
  auto task = mode == Dissimilarity::SCORE ? EDLIB_TASK_DISTANCE : EDLIB_TASK_PATH;
  auto config = edlibNewAlignConfig(-1, EDLIB_MODE_NW, task, NULL, 0);
  return distance_edlib(a, b, config, mode);
}

double align_edlib_extend(
    const std::string a,
    const std::string b,
    std::string dissimilarity
) {
  Dissimilarity mode = parse_dissimilarity(dissimilarity);
  auto task = mode == Dissimilarity::SCORE ? EDLIB_TASK_DISTANCE : EDLIB_TASK_PATH;
  auto config = edlibNewAlignConfig(-1, EDLIB_MODE_SHW, task, NULL, 0);
  // SHW is (query, target) with the query's suffix free; use the shorter
  // sequence as query so all three backends share L = min(l1, l2).
  if (a.size() <= b.size()) {
    return distance_edlib(a, b, config, mode);
  }
  return distance_edlib(b, a, config, mode);
}

template<enum AlignmentSpan span>
double distance_ksw2(
    const std::string &a,
    const std::string &b,
    Ksw2Aligner &aligner,
    Dissimilarity dissimilarity
) {
  if (!aligner.align(a, b, span)) return 1.0;
  if (dissimilarity == Dissimilarity::SCORE) {
    int L = score_ref_length(a.size(), b.size(), span);
    return ksw_score_distance(
        aligner.getAlignmentScore(), L, aligner.match, aligner.mismatch
    );
  }
  std::string cigar = aligner.getCIGAR();
  if (cigar.empty()) return 1.0;
  if constexpr (span == AlignmentSpan::EXTEND) {
    return distance_from_cigar_extend(cigar);
  } else {
    return distance_from_cigar(cigar);
  }
}

template double distance_ksw2<AlignmentSpan::GLOBAL>(
    const std::string &a, const std::string &b, Ksw2Aligner &aligner,
    Dissimilarity dissimilarity);
template double distance_ksw2<AlignmentSpan::EXTEND>(
    const std::string &a, const std::string &b, Ksw2Aligner &aligner,
    Dissimilarity dissimilarity);

template<enum AlignmentSpan span>
double distance_ksw2(
    const SequenceView &a,
    const SequenceView &b,
    Ksw2Aligner &aligner,
    Dissimilarity dissimilarity
) {
  if (!aligner.align(a, b, span)) return 1.0;
  if (dissimilarity == Dissimilarity::SCORE) {
    int L = score_ref_length(a.size(), b.size(), span);
    return ksw_score_distance(
        aligner.getAlignmentScore(), L, aligner.match, aligner.mismatch
    );
  }
  std::string cigar = aligner.getCIGAR();
  if (cigar.empty()) return 1.0;
  if constexpr (span == AlignmentSpan::EXTEND) {
    return distance_from_cigar_extend(cigar);
  } else {
    return distance_from_cigar(cigar);
  }
}

template double distance_ksw2<AlignmentSpan::GLOBAL>(
    const SequenceView &a, const SequenceView &b, Ksw2Aligner &aligner,
    Dissimilarity dissimilarity);
template double distance_ksw2<AlignmentSpan::EXTEND>(
    const SequenceView &a, const SequenceView &b, Ksw2Aligner &aligner,
    Dissimilarity dissimilarity);

template<enum AlignmentSpan span>
std::string cigar_ksw2(const std::string &a, const std::string &b,
                       Ksw2Aligner &aligner) {
  if (!aligner.align(a, b, span)) return "";
  return aligner.getCIGAR();
}

template std::string cigar_ksw2<AlignmentSpan::GLOBAL>(
    const std::string &a, const std::string &b, Ksw2Aligner &aligner);
template std::string cigar_ksw2<AlignmentSpan::EXTEND>(
    const std::string &a, const std::string &b, Ksw2Aligner &aligner);

std::string cigar_ksw2_global(const std::string &a, const std::string &b,
                              int match, int mismatch,
                              int gap_open, int gap_extend,
                              int gap_open2, int gap_extend2) {
  Ksw2Aligner aligner{match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2};
  return cigar_ksw2<AlignmentSpan::GLOBAL>(a, b, aligner);
}

std::string cigar_ksw2_extend(const std::string &a, const std::string &b,
                              int match, int mismatch,
                              int gap_open, int gap_extend,
                              int gap_open2, int gap_extend2) {
  Ksw2Aligner aligner{match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2};
  return cigar_ksw2<AlignmentSpan::EXTEND>(a, b, aligner);
}

template<enum AlignmentSpan span>
std::pair<double, std::string> distance_and_cigar_ksw2(
    const std::string &a,
    const std::string &b,
    Ksw2Aligner &aligner,
    Dissimilarity dissimilarity
) {
  if (!aligner.align(a, b, span)) return {1.0, ""};
  std::string cigar = aligner.getCIGAR();
  if (cigar.empty()) return {1.0, ""};
  double d;
  if (dissimilarity == Dissimilarity::SCORE) {
    int L = score_ref_length(a.size(), b.size(), span);
    d = ksw_score_distance(
        aligner.getAlignmentScore(), L, aligner.match, aligner.mismatch
    );
  } else if constexpr (span == AlignmentSpan::EXTEND) {
    d = distance_from_cigar_extend(cigar);
  } else {
    d = distance_from_cigar(cigar);
  }
  return {d, cigar};
}

template std::pair<double, std::string> distance_and_cigar_ksw2<AlignmentSpan::GLOBAL>(
    const std::string &a, const std::string &b, Ksw2Aligner &aligner,
    Dissimilarity dissimilarity);
template std::pair<double, std::string> distance_and_cigar_ksw2<AlignmentSpan::EXTEND>(
    const std::string &a, const std::string &b, Ksw2Aligner &aligner,
    Dissimilarity dissimilarity);

double align_ksw2_global(const std::string a, const std::string b,
             int match, int mismatch,
             int gap_open, int gap_extend,
             int gap_open2, int gap_extend2,
             std::string dissimilarity) {
  Dissimilarity mode = parse_dissimilarity(dissimilarity);
  Ksw2Aligner aligner{match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2};
  aligner.setScoreOnly(mode == Dissimilarity::SCORE);
  return distance_ksw2<AlignmentSpan::GLOBAL>(a, b, aligner, mode);
}

double align_ksw2_extend(const std::string a, const std::string b,
             int match, int mismatch,
             int gap_open, int gap_extend,
             int gap_open2, int gap_extend2,
             std::string dissimilarity) {
  Dissimilarity mode = parse_dissimilarity(dissimilarity);
  Ksw2Aligner aligner{match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2};
  aligner.setScoreOnly(mode == Dissimilarity::SCORE);
  return distance_ksw2<AlignmentSpan::EXTEND>(a, b, aligner, mode);
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
    // dist_wfa2() match is a non-negative bonus; WFA2 requires <= 0.
    match = -match;
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
  } else if (dist_method == "ksw2") {
    int match = element_as_int(dist_config, "match", "dist_config");
    int mismatch = element_as_int(dist_config, "mismatch", "dist_config");
    int gap_open = element_as_int(dist_config, "gap_open", "dist_config");
    int gap_extend = element_as_int(dist_config, "gap_extend", "dist_config");
    int gap_open2 = element_as_int(dist_config, "gap_open2", "dist_config");
    int gap_extend2 = element_as_int(dist_config, "gap_extend2", "dist_config");
    if (span == 0) {
      return align_from_compressed_cigar(
        a, b,
        cigar_ksw2_global(
          a, b, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2
        )
      );
    } else if (span == 1) {
      return align_from_compressed_cigar(
        a, b,
        cigar_ksw2_extend(
          a, b, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2
        )
      );
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

#ifdef OPTIMOTU_R
Rcpp::RObject add_gapstats(Rcpp::DataFrame df, std::string cigar_column) {
  Rcpp::CharacterVector cigars = df[cigar_column];
  std::vector<std::string> cigars_vec = Rcpp::as<std::vector<std::string>>(cigars);
  std::vector<int> align_length;
  std::vector<int> n_insert;
  std::vector<int> n_delete;
  std::vector<int> max_insert;
  std::vector<int> max_delete;
  for (const std::string & cigar : cigars_vec) {
    auto [dist, alen, n_ins, n_del, max_ins, max_del] =
      dist_gapstats_from_cigar(cigar);
    align_length.push_back(alen);
    n_insert.push_back(n_ins);
    n_delete.push_back(n_del);
    max_insert.push_back(max_ins);
    max_delete.push_back(max_del);
  }
  df["align_length"] = align_length;
  df["n_insert"] = n_insert;
  df["n_delete"] = n_delete;
  df["max_insert"] = max_insert;
  df["max_delete"] = max_delete;
  df.attr("class") = Rcpp::CharacterVector::create("tbl_df", "tbl", "data.frame");
  df.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, -align_length.size());
  return df;
}
#endif //OPTIMOTU_R
