// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_HAMMING_STRIPE_H
#define OPTIMOTU_HAMMING_STRIPE_H

#include "hamming_kernel.h"
#include <cstdint>

// Score one packed query against a word-major stripe of targets
// (ignore_gaps / pdistB semantics: invalid bases are already zero).
void hamming_score_stripe(
  const uint64_t *query_seq,
  const uint64_t *query_mask,
  const uint64_t *seq_major,
  const uint64_t *mask_major,
  int ulen,
  int mulen,
  int n_target,
  int n_cols,
  int min_overlap,
  double *dist_out,
  int *n_ok_out,
  int *n_match_out
);

// Pairwise pdistB over all packed words (terminal-gap zeros are free).
int hamming_pdist_pair(
  const uint64_t *a,
  const uint64_t *ma,
  const uint64_t *b,
  const uint64_t *mb,
  int ulen,
  int mulen,
  int min_overlap,
  int *num_ok,
  int *num_matches,
  double *dist
);

#endif
