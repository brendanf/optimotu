#ifndef OPTIMOTU_PACKEDSEQUENCESET_H
#define OPTIMOTU_PACKEDSEQUENCESET_H

#include <vector>
#include <string>
#include <cstdint>
#include <utility>
#include <tuple>
#include "optimotu.h"
#include "SequenceView.h"

// Analogous to SequenceSetB in C code
struct PackedSequenceSet {
  int num_seqs = 0, alen = 0, ulen = 0, mulen = 0;
  std::vector<std::vector<uint64_t>> packed_seq;
  std::vector<std::vector<uint64_t>> mask;
  std::vector<int> start;
  std::vector<int> end;
  PackedSequenceSet(const SequenceSet &seq);
  PackedSequenceSet(
    const SequenceSet & seq1,
    const SequenceSet & seq2
  );
  // Convenience for callers that still own std::string buffers.
  PackedSequenceSet(const std::vector<std::string> &seq);
  PackedSequenceSet(
    const std::vector<std::string> & seq1,
    const std::vector<std::string> & seq2
  );

  // Approximate heap bytes for packed bits + mask + start/end.
  static std::size_t estimate_bytes(std::size_t n, int alen);

  double dist(
    const int i,
    const int j,
    const int min_overlap,
    const bool ignore_gap
  ) const;
  std::pair<int, double> score_and_dist(
    const int i,
    const int j,
    const int min_overlap,
    const bool ignore_gap
  ) const;
  std::tuple<bool, int, double> success_score_and_dist(
    const int i,
    const int j,
    const int min_overlap,
    const bool ignore_gap
  ) const;
  std::tuple<bool, int, double, int, int, int, int> success_score_and_dist_gap(
    const int i,
    const int j,
    const int min_overlap,
    const bool ignore_gap
  ) const;
  bool verify(const SequenceSet & seq) const;

  void fill_word_major(
    int begin,
    int end,
    std::vector<uint64_t> &seq_out,
    std::vector<uint64_t> &mask_out
  ) const;
  void fill_word_major_idx(
    const int *idx,
    int n_idx,
    std::vector<uint64_t> &seq_out,
    std::vector<uint64_t> &mask_out
  ) const;

  // ignore_gaps / pdistB only. seq_major and mask_major from fill_word_major*.
  void score_query_vs_word_major(
    int query_i,
    const uint64_t *seq_major,
    const uint64_t *mask_major,
    int n_target,
    int n_cols,
    int min_overlap,
    double *dist_out,
    int *n_ok_out,
    int *n_match_out = nullptr
  ) const;
};

#endif
