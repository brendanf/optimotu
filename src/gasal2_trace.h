// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_GASAL2_TRACE_H
#define OPTIMOTU_GASAL2_TRACE_H

#include "Ksw2Aligner.h"
#include "gasal2_edit.h"
#include "pairwise_alignment.h"

#include <edlib.h>

struct Gasal2HostAlign {
  bool edit;
  int match;
  int mismatch;
  int gap_open;
  int gap_extend;
  Ksw2Aligner ksw;
  EdlibAlignConfig edlib;

  Gasal2HostAlign(
    int match_,
    int mismatch_,
    int gap_open_,
    int gap_extend_
  ) : edit(gasal2_uses_edit(match_, mismatch_, gap_open_, gap_extend_)),
      match(match_),
      mismatch(mismatch_),
      gap_open(gap_open_),
      gap_extend(gap_extend_),
      ksw(match_, mismatch_, gap_open_, gap_extend_, gap_open_, gap_extend_),
      edlib(edlibNewAlignConfig(
        -1,
        EDLIB_MODE_NW,
        EDLIB_TASK_PATH,
        nullptr,
        0
      )) {}

  bool pair_ok(
    Dissimilarity dissimilarity,
    double l1,
    double l2,
    double threshold
  ) const {
    if (edit) {
      return pair_length_ok(
        dissimilarity,
        AlignmentSpan::GLOBAL,
        l1,
        l2,
        threshold,
        0,
        mismatch,
        gap_open,
        gap_extend,
        gap_open,
        gap_extend
      );
    }
    return pair_length_ok(
      dissimilarity,
      AlignmentSpan::GLOBAL,
      l1,
      l2,
      threshold,
      0,
      ksw_effective_match(match) + ksw_effective_mismatch(mismatch),
      gap_open,
      gap_extend,
      gap_open,
      gap_extend
    );
  }

  template<typename Seq>
  double identity(const Seq &query, const Seq &target) {
    if (edit) {
      return distance_edlib(
        query, target, edlib, Dissimilarity::IDENTITY
      );
    }
    ksw.setScoreOnly(false);
    return distance_ksw2(query, target, ksw, Dissimilarity::IDENTITY);
  }
};

#endif
