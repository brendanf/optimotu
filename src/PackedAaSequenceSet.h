// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_PACKEDAASEQUENCESET_H
#define OPTIMOTU_PACKEDAASEQUENCESET_H

#include "SequenceView.h"
#include <cstddef>
#include <cstdint>
#include <tuple>
#include <vector>

// Packed equal-length amino-acid alignments.
// Code 0 = gap; 1-26 = A-Z; 27 = *; 28 = unknown. Lowercase is skipped
// as an A2M insertion (same rule as PackedSequenceSet).
struct PackedAaSequenceSet {
  static constexpr int kLut = 32;

  int num_seqs = 0;
  int alen = 0;
  std::vector<uint8_t> codes;
  std::vector<int> start;
  std::vector<int> end;

  PackedAaSequenceSet() = default;
  explicit PackedAaSequenceSet(const SequenceSet &seq);
  PackedAaSequenceSet(const SequenceSet &seq1, const SequenceSet &seq2);
  explicit PackedAaSequenceSet(const std::vector<std::string> &seq);
  PackedAaSequenceSet(
    const std::vector<std::string> &seq1,
    const std::vector<std::string> &seq2
  );

  static std::size_t estimate_bytes(std::size_t n, int alen);

  const uint8_t *row(int i) const {
    return codes.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(alen);
  }

  void fill_pos_major(int begin, int end, std::vector<uint8_t> &out) const;
  void fill_pos_major_idx(
    const int *idx,
    int n_idx,
    std::vector<uint8_t> &out
  ) const;

  double dist(
    int i,
    int j,
    int min_overlap,
    bool ignore_gaps,
    const float *weight_lut
  ) const;

  std::tuple<bool, int, double> success_score_and_dist(
    int i,
    int j,
    int min_overlap,
    bool ignore_gaps,
    const float *weight_lut
  ) const;

  void score_query_vs_pos_major(
    int query_i,
    const uint8_t *pos_major,
    const int *target_idx,
    int n_target,
    int n_cols,
    int min_overlap,
    bool ignore_gaps,
    const float *weight_lut,
    double dist_threshold,
    bool early_exit,
    double *dist_out,
    int *n_ok_out
  ) const;
};

uint8_t pack_aa_char(char c);

#endif
