#ifndef OPTIMOTU_WFA_IDENTITY_BOUND_H
#define OPTIMOTU_WFA_IDENTITY_BOUND_H

#include <algorithm>
#include <cmath>
#include <limits>

// Band and exclusive WFA score cap so that no global alignment with CIGAR
// identity distance <= threshold can fall outside the band or exceed the
// cap. Derivation: docs/wfa-identity-score-bound.md

struct WfaIdentityBound {
  int min_k = 0;
  int max_k = 0;
  // Pass to WFAligner::setMaxAlignmentSteps (abort when score >= this).
  int max_alignment_steps = 1;
};

// Max affine/dual-affine cost of G gap characters (worst-case runs).
inline int wfa_max_gap_cost(
    int n_gap_chars,
    int gap_open,
    int gap_extend,
    int gap_open2,
    int gap_extend2
) {
  if (n_gap_chars <= 0) return 0;
  long long g = n_gap_chars;
  long long unit = std::max(
      (long long)gap_open + gap_extend,
      (long long)gap_open2 + gap_extend2
  );
  long long many = unit * g;
  long long one = std::max(
      (long long)gap_open + (long long)gap_extend * g,
      (long long)gap_open2 + (long long)gap_extend2 * g
  );
  long long c = std::max(many, one);
  if (c > std::numeric_limits<int>::max()) {
    return std::numeric_limits<int>::max();
  }
  return (int)c;
}

inline double wfa_identity_max_edits(
    double l1,
    double l2,
    double threshold
) {
  if (threshold <= 0.0) return 0.0;
  return threshold * (l1 + l2) / (2.0 - threshold);
}

// Inclusive max WFA score of any identity-feasible CIGAR (match<=0 ignored).
inline int wfa_identity_max_score(
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
  double delta = l2 - l1;
  if (delta < 0.0) delta = 0.0;
  int g_delta = (int)std::ceil(delta);

  long long s;
  if (threshold <= 0.0) {
    s = wfa_max_gap_cost(
        g_delta, gap_open, gap_extend, gap_open2, gap_extend2
    );
  } else {
    double x_max = threshold * l2 - delta;
    if (x_max < 0.0) x_max = 0.0;
    int X = (int)std::ceil(x_max);
    int g_max = (int)std::ceil(
        wfa_identity_max_edits(l1, l2, threshold)
    );
    long long s1 = (long long)mismatch * X +
      wfa_max_gap_cost(
          g_delta, gap_open, gap_extend, gap_open2, gap_extend2
      );
    long long s2 = wfa_max_gap_cost(
        g_max, gap_open, gap_extend, gap_open2, gap_extend2
    );
    s = std::max(s1, s2);
  }
  // Positive match is a per-column cost; match<=0 cannot raise the cap.
  if (match > 0) {
    long long m = (long long)std::ceil(l1);
    if (m < 0) m = 0;
    s += (long long)match * m;
  }
  if (s > std::numeric_limits<int>::max()) {
    return std::numeric_limits<int>::max();
  }
  if (s < 0) return 0;
  return (int)s;
}

inline WfaIdentityBound wfa_identity_bound(
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
  WfaIdentityBound out;
  double sigma = 1.0 - threshold;
  double denom = 1.0 + sigma;
  if (denom <= 0.0) denom = 1.0;
  out.max_k = (int)std::ceil((l2 - l1 * sigma) / denom);
  if (extend_span) {
    out.min_k = -(int)std::ceil(l1 * sigma);
  } else {
    out.min_k = -(int)std::ceil((l1 - l2 * sigma) / denom);
  }
  int smax = wfa_identity_max_score(
      l1, l2, threshold, match, mismatch,
      gap_open, gap_extend, gap_open2, gap_extend2
  );
  if (smax >= std::numeric_limits<int>::max()) {
    out.max_alignment_steps = std::numeric_limits<int>::max();
  } else {
    out.max_alignment_steps = smax + 1;
  }
  if (out.max_alignment_steps < 1) out.max_alignment_steps = 1;
  return out;
}

#endif
