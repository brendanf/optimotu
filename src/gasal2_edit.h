// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_GASAL2_EDIT_H
#define OPTIMOTU_GASAL2_EDIT_H

#include "score_dissimilarity.h"

// Unit (or uniformly scaled) Levenshtein: zero match bonus, linear
// gaps, and equal substitution/indel costs. Affine or unequal linear
// costs stay on the KSW2 similarity path (match 0 maps to +1).
inline bool gasal2_uses_edit(
  int match,
  int mismatch,
  int gap_open,
  int gap_extend
) {
  return match <= 0 && gap_open == 0 && mismatch > 0 && gap_extend > 0 &&
    mismatch == gap_extend;
}

inline int gasal2_gpu_match(
  int match,
  int mismatch,
  int gap_open,
  int gap_extend
) {
  if (gasal2_uses_edit(match, mismatch, gap_open, gap_extend)) {
    return 0;
  }
  return match > 0 ? match : 1;
}

inline int gasal2_gpu_mismatch(int mismatch) {
  return mismatch > 0 ? -mismatch : -1;
}

inline double gasal2_score_distance(
  int score,
  int L,
  int match,
  int mismatch,
  int gap_open,
  int gap_extend
) {
  if (gasal2_uses_edit(match, mismatch, gap_open, gap_extend)) {
    // Kernel returns similarity with match 0, so -score is the edit
    // (or uniformly scaled edit) cost.
    return cost_score_distance(-score, L, 0, mismatch);
  }
  return ksw_score_distance(score, L, match, mismatch);
}

#endif
