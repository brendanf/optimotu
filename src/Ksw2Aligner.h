// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_KSW2ALIGNER_H
#define OPTIMOTU_KSW2ALIGNER_H

#include "alignment_enums.h"
#include "SequenceView.h"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

extern "C" {
#include "ksw2.h"
}

// Thin wrapper around KSW2 global / extension alignment with affine or
// dual-affine gaps. Match/mismatch/gap parameters use the same convention as
// dist_wfa2() / dist_ksw2() in R: positive match is a bonus, positive
// mismatch and gap values are penalties.
class Ksw2Aligner {
public:
  Ksw2Aligner(
    int match = 0,
    int mismatch = 1,
    int gap_open = 0,
    int gap_extend = 1,
    int gap_open2 = 0,
    int gap_extend2 = 1
  );

  ~Ksw2Aligner();

  Ksw2Aligner(const Ksw2Aligner &) = delete;
  Ksw2Aligner &operator=(const Ksw2Aligner &) = delete;

  // Symmetric band width. Use w < 0 to disable banding.
  void setBandWidth(int w) { band_width = w; }

  // Align a (query) vs b (target). Callers should pass the longer sequence as
  // query and the shorter as target when following the WFA2 worker convention.
  // Returns true on a completed path; false if z-dropped / incomplete.
  bool align(
    const char *a,
    int a_len,
    const char *b,
    int b_len,
    AlignmentSpan span = AlignmentSpan::GLOBAL
  );

  bool align(
    const std::string &a,
    const std::string &b,
    AlignmentSpan span = AlignmentSpan::GLOBAL
  ) {
    return align(
      a.data(), static_cast<int>(a.size()),
      b.data(), static_cast<int>(b.size()),
      span
    );
  }

  bool align(
    const SequenceView &a,
    const SequenceView &b,
    AlignmentSpan span = AlignmentSpan::GLOBAL
  ) {
    return align(
      a.data(), static_cast<int>(a.size()),
      b.data(), static_cast<int>(b.size()),
      span
    );
  }

  // Compressed CIGAR using =/X/I/D. Empty if the last align() failed.
  std::string getCIGAR() const;

  int getAlignmentScore() const { return last_score; }

  int match = 0;
  int mismatch = 1;
  int gap_open = 0;
  int gap_extend = 1;
  int gap_open2 = 0;
  int gap_extend2 = 1;

private:
  void rebuild_matrix();
  void encode(const char *seq, int len, std::vector<uint8_t> &out) const;
  void free_cigar();

  int8_t mat[25]{};
  int8_t sc_match = 1;
  int8_t sc_mismatch = -1;
  bool dual_affine = false;
  int band_width = -1;
  int last_score = 0;
  bool last_ok = false;

  ksw_extz_t ez{};
  std::vector<uint8_t> query_enc;
  std::vector<uint8_t> target_enc;
  // Encoded sequences retained for EQ/X CIGAR conversion
  std::vector<uint8_t> last_query;
  std::vector<uint8_t> last_target;
};

#endif
