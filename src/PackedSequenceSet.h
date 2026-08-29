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
  int num_seqs, alen, ulen, mulen;
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
};

#endif
