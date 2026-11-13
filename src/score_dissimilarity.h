#ifndef OPTIMOTU_SCORE_DISSIMILARITY_H
#define OPTIMOTU_SCORE_DISSIMILARITY_H

#include "alignment_enums.h"
#include "wfa_identity_bound.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

inline int score_ref_length(double l1, double l2, AlignmentSpan span) {
  if (l1 < 0) l1 = 0;
  if (l2 < 0) l2 = 0;
  if (span == AlignmentSpan::EXTEND) {
    return static_cast<int>(std::min(l1, l2));
  }
  return static_cast<int>(std::max(l1, l2));
}

inline int edlib_score_k(double L, double threshold) {
  int k = static_cast<int>(std::floor(threshold * L));
  if (k < 0) k = 0;
  return k;
}

// Cost-minimizing aligners (WFA2, Edlib). match is the match cost the
// aligner received (0, or negative if a bonus was mapped to a cost).
// WFA reports a Smith-Waterman-style classic score. With match cost 0
// that is the negation of the internal cost; undo it before the
// substitution-unit formula.
inline int wfa_reported_to_cost(int reported, int match_cost) {
  if (match_cost <= 0 && reported < 0) return -reported;
  return reported;
}

inline double cost_score_distance(
    int score,
    int L,
    int match_cost,
    int mismatch
) {
  if (L <= 0) return 1.0;
  int delta = mismatch - match_cost;
  if (delta <= 0) return 1.0;
  double d = (static_cast<double>(score) -
              static_cast<double>(match_cost) * static_cast<double>(L)) /
             (static_cast<double>(delta) * static_cast<double>(L));
  if (d < 0.0) d = 0.0;
  if (d > 1.0) d = 1.0;
  if (!std::isfinite(d)) return 1.0;
  return d;
}

// KSW2 similarity. match==0 maps to mu=1 as in Ksw2Aligner.
inline int ksw_effective_match(int match) {
  return match > 0 ? match : 1;
}

inline int ksw_effective_mismatch(int mismatch) {
  return mismatch > 0 ? mismatch : 1;
}

inline double ksw_score_distance(int similarity, int L, int match, int mismatch) {
  if (L <= 0) return 1.0;
  int mu = ksw_effective_match(match);
  int nu = ksw_effective_mismatch(mismatch);
  int delta = mu + nu;
  if (delta <= 0) return 1.0;
  double d = (static_cast<double>(mu) * static_cast<double>(L) -
              static_cast<double>(similarity)) /
             (static_cast<double>(delta) * static_cast<double>(L));
  if (d < 0.0) d = 0.0;
  if (d > 1.0) d = 1.0;
  if (!std::isfinite(d)) return 1.0;
  return d;
}

inline WfaIdentityBound ksw_score_bound(
    double l1,
    double l2,
    double threshold,
    int match,
    int mismatch,
    int gap_open,
    int gap_extend,
    int gap_open2,
    int gap_extend2,
    bool extend_span = false
) {
  int mu = ksw_effective_match(match);
  int nu = ksw_effective_mismatch(mismatch);
  return wfa_score_bound(
      l1, l2, threshold, 0, mu + nu,
      gap_open, gap_extend, gap_open2, gap_extend2, extend_span
  );
}

inline WfaIdentityBound choose_alignment_bound(
    Dissimilarity mode,
    double l1,
    double l2,
    double threshold,
    int match,
    int mismatch,
    int gap_open,
    int gap_extend,
    int gap_open2,
    int gap_extend2,
    bool extend_span = false,
    bool ksw_similarity = false
) {
  if (mode == Dissimilarity::SCORE) {
    if (ksw_similarity) {
      return ksw_score_bound(
          l1, l2, threshold, match, mismatch,
          gap_open, gap_extend, gap_open2, gap_extend2, extend_span
      );
    }
    return wfa_score_bound(
        l1, l2, threshold, match, mismatch,
        gap_open, gap_extend, gap_open2, gap_extend2, extend_span
    );
  }
  return wfa_identity_bound(
      l1, l2, threshold, match, mismatch,
      gap_open, gap_extend, gap_open2, gap_extend2, extend_span
  );
}

inline bool pair_length_ok(
    Dissimilarity mode,
    AlignmentSpan span,
    double l1,
    double l2,
    double threshold,
    int match,
    int mismatch,
    int gap_open,
    int gap_extend,
    int gap_open2,
    int gap_extend2
) {
  if (l2 < l1) {
    double tmp = l1;
    l1 = l2;
    l2 = tmp;
  }
  if (l2 <= 0.0) return false;
  if (mode == Dissimilarity::SCORE) {
    if (span == AlignmentSpan::EXTEND) return true;
    return score_global_length_ok(
        l1, l2, threshold, match, mismatch,
        gap_open, gap_extend, gap_open2, gap_extend2
    );
  }
  return (l1 / l2) >= (1.0 - threshold);
}

#endif
