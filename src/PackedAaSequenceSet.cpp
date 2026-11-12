// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "PackedAaSequenceSet.h"
#include "optimotu.h"
#include <algorithm>
#include <cstring>
#include <vector>

#if defined(__AVX2__)
#include <immintrin.h>
#elif defined(__SSE2__)
#include <emmintrin.h>
#endif

uint8_t pack_aa_char(char c) {
  if (c >= 'a' && c <= 'z') {
    return 255;
  }
  if (c == '-' || c == '.' || c == ' ') {
    return 0;
  }
  if (c >= 'A' && c <= 'Z') {
    return static_cast<uint8_t>(c - 'A' + 1);
  }
  if (c == '*') {
    return 27;
  }
  return 28;
}

namespace {

int packed_aa_width(const SequenceView &s) {
  int n = 0;
  for (std::size_t k = 0; k < s.size(); ++k) {
    if (s[k] < 'a' || s[k] > 'z') {
      ++n;
    }
  }
  return n;
}

void pack_one(const SequenceView &s, uint8_t *out, int alen, int *start, int *end) {
  int n = 0;
  for (std::size_t k = 0; k < s.size(); ++k) {
    const uint8_t code = pack_aa_char(s[k]);
    if (code == 255) {
      continue;
    }
    if (n >= alen) {
      OPTIMOTU_STOP(
        "Amino-acid sequences must have equal length after removing "
        "lowercase insertions"
      );
    }
    out[n++] = code;
  }
  if (n != alen) {
    OPTIMOTU_STOP(
      "Amino-acid sequences must have equal length after removing "
      "lowercase insertions"
    );
  }
  int s0 = 0;
  while (s0 < alen && out[s0] == 0) {
    ++s0;
  }
  int e0 = alen;
  while (e0 > s0 && out[e0 - 1] == 0) {
    --e0;
  }
  *start = s0;
  *end = e0;
}

void pack_set(
  const SequenceSet &seq,
  std::vector<uint8_t> &codes,
  std::vector<int> &start,
  std::vector<int> &end,
  int *num_seqs,
  int *alen
) {
  *num_seqs = static_cast<int>(seq.size());
  if (*num_seqs == 0) {
    *alen = 0;
    codes.clear();
    start.clear();
    end.clear();
    return;
  }
  *alen = packed_aa_width(seq[0]);
  codes.assign(
    static_cast<std::size_t>(*num_seqs) * static_cast<std::size_t>(*alen),
    0
  );
  start.assign(*num_seqs, 0);
  end.assign(*num_seqs, 0);
  for (int i = 0; i < *num_seqs; ++i) {
    pack_one(
      seq[static_cast<std::size_t>(i)],
      codes.data() + static_cast<std::size_t>(i) * static_cast<std::size_t>(*alen),
      *alen,
      &start[static_cast<std::size_t>(i)],
      &end[static_cast<std::size_t>(i)]
    );
  }
}

inline void identity_site(
  uint8_t a,
  uint8_t b,
  bool ignore_gaps,
  int *n_ok,
  int *n_match
) {
  if (a == 0 && b == 0) {
    return;
  }
  if (ignore_gaps && (a == 0 || b == 0)) {
    return;
  }
  ++(*n_ok);
  if (a == b) {
    ++(*n_match);
  }
}

inline void identity_site_u16(
  uint8_t a,
  uint8_t b,
  bool ignore_gaps,
  uint16_t *n_ok,
  uint16_t *n_match
) {
  int ok = 0;
  int match = 0;
  identity_site(a, b, ignore_gaps, &ok, &match);
  *n_ok = static_cast<uint16_t>(*n_ok + ok);
  *n_match = static_cast<uint16_t>(*n_match + match);
}

#if defined(__SSE2__)
inline __m128i loadu_si128(const uint8_t *p) {
  return _mm_loadu_si128(reinterpret_cast<const __m128i *>(p));
}

inline __m128i skip_mask_16(uint8_t a, __m128i vb, bool ignore_gaps) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i b_gap = _mm_cmpeq_epi8(vb, zero);
  if (ignore_gaps) {
    if (a == 0) {
      return _mm_set1_epi8(static_cast<char>(0xFF));
    }
    return b_gap;
  }
  if (a == 0) {
    return b_gap;
  }
  return zero;
}

inline void add_u16x16(uint16_t *dst, __m128i bytes_0_or_1) {
  const __m128i zero = _mm_setzero_si128();
  const __m128i lo = _mm_unpacklo_epi8(bytes_0_or_1, zero);
  const __m128i hi = _mm_unpackhi_epi8(bytes_0_or_1, zero);
  const __m128i dlo = _mm_loadu_si128(reinterpret_cast<const __m128i *>(dst));
  const __m128i dhi = _mm_loadu_si128(
    reinterpret_cast<const __m128i *>(dst + 8)
  );
  _mm_storeu_si128(
    reinterpret_cast<__m128i *>(dst),
    _mm_add_epi16(dlo, lo)
  );
  _mm_storeu_si128(
    reinterpret_cast<__m128i *>(dst + 8),
    _mm_add_epi16(dhi, hi)
  );
}

inline void identity_column_16(
  uint8_t a,
  const uint8_t *col,
  bool ignore_gaps,
  uint16_t *n_ok,
  uint16_t *n_match
) {
  const __m128i va = _mm_set1_epi8(static_cast<char>(a));
  const __m128i vb = loadu_si128(col);
  const __m128i skip = skip_mask_16(a, vb, ignore_gaps);
  const __m128i ones = _mm_set1_epi8(static_cast<char>(0xFF));
  const __m128i one = _mm_set1_epi8(1);
  const __m128i valid = _mm_andnot_si128(skip, ones);
  const __m128i eq = _mm_and_si128(_mm_cmpeq_epi8(va, vb), valid);
  add_u16x16(n_ok, _mm_and_si128(valid, one));
  add_u16x16(n_match, _mm_and_si128(eq, one));
}

inline void identity_pair_16(
  const uint8_t *a,
  const uint8_t *b,
  bool ignore_gaps,
  int *n_ok,
  int *n_match
) {
  const __m128i va = loadu_si128(a);
  const __m128i vb = loadu_si128(b);
  const __m128i zero = _mm_setzero_si128();
  const __m128i a_gap = _mm_cmpeq_epi8(va, zero);
  const __m128i b_gap = _mm_cmpeq_epi8(vb, zero);
  const __m128i skip = ignore_gaps
    ? _mm_or_si128(a_gap, b_gap)
    : _mm_and_si128(a_gap, b_gap);
  const __m128i ones = _mm_set1_epi8(static_cast<char>(0xFF));
  const __m128i valid = _mm_andnot_si128(skip, ones);
  const __m128i eq = _mm_and_si128(_mm_cmpeq_epi8(va, vb), valid);
  *n_ok += __builtin_popcount(_mm_movemask_epi8(valid));
  *n_match += __builtin_popcount(_mm_movemask_epi8(eq));
}
#endif

#if defined(__AVX2__)
inline __m128i loadu_si64(const uint8_t *p) {
  long long tmp;
  std::memcpy(&tmp, p, sizeof(tmp));
  return _mm_cvtsi64_si128(tmp);
}

inline void weighted_column_8(
  const float *prof_row,
  uint8_t a,
  const uint8_t *col,
  bool ignore_gaps,
  int *n_ok,
  double *acc
) {
  const __m128i codes = loadu_si64(col);
  const __m128i zero8 = _mm_setzero_si128();
  const __m128i b_gap = _mm_cmpeq_epi8(codes, zero8);
  __m128i skip;
  if (ignore_gaps) {
    skip = a == 0 ? _mm_set1_epi8(static_cast<char>(0xFF)) : b_gap;
  } else if (a == 0) {
    skip = b_gap;
  } else {
    skip = zero8;
  }
  const __m128i ones = _mm_set1_epi8(static_cast<char>(0xFF));
  // High 8 bytes of the 64-bit load are zero; keep only the low 8 lanes.
  const __m128i low8 = _mm_cvtsi64_si128(-1LL);
  const __m128i valid8 = _mm_and_si128(_mm_andnot_si128(skip, ones), low8);
  const __m256i valid32 = _mm256_cvtepu8_epi32(valid8);
  const __m256i is_valid = _mm256_cmpgt_epi32(
    valid32,
    _mm256_setzero_si256()
  );
  const __m256i one = _mm256_and_si256(is_valid, _mm256_set1_epi32(1));
  const __m256i nok = _mm256_add_epi32(
    _mm256_loadu_si256(reinterpret_cast<const __m256i *>(n_ok)),
    one
  );
  _mm256_storeu_si256(reinterpret_cast<__m256i *>(n_ok), nok);

  const __m256i idx = _mm256_cvtepu8_epi32(codes);
  const __m256 gathered = _mm256_i32gather_ps(prof_row, idx, 4);
  const __m256 masked = _mm256_and_ps(
    gathered,
    _mm256_castsi256_ps(is_valid)
  );
  const __m128 wlo = _mm256_castps256_ps128(masked);
  const __m128 whi = _mm256_extractf128_ps(masked, 1);
  _mm256_storeu_pd(
    acc,
    _mm256_add_pd(_mm256_loadu_pd(acc), _mm256_cvtps_pd(wlo))
  );
  _mm256_storeu_pd(
    acc + 4,
    _mm256_add_pd(_mm256_loadu_pd(acc + 4), _mm256_cvtps_pd(whi))
  );
}
#endif

inline void weighted_site(
  uint8_t a,
  uint8_t b,
  bool ignore_gaps,
  const float *lut,
  int *n_ok,
  double *acc
) {
  if (a == 0 && b == 0) {
    return;
  }
  if (ignore_gaps && (a == 0 || b == 0)) {
    return;
  }
  ++(*n_ok);
  *acc += static_cast<double>(
    lut[static_cast<std::size_t>(a) * PackedAaSequenceSet::kLut + b]
  );
}

}  // namespace

std::size_t PackedAaSequenceSet::estimate_bytes(std::size_t n, int alen) {
  if (n == 0 || alen <= 0) {
    return 0;
  }
  return n * (static_cast<std::size_t>(alen) + 2 * sizeof(int));
}

PackedAaSequenceSet::PackedAaSequenceSet(const SequenceSet &seq) {
  pack_set(seq, codes, start, end, &num_seqs, &alen);
}

PackedAaSequenceSet::PackedAaSequenceSet(
  const SequenceSet &seq1,
  const SequenceSet &seq2
) {
  SequenceSet both;
  both.reserve(seq1.size() + seq2.size());
  both.insert(both.end(), seq1.begin(), seq1.end());
  both.insert(both.end(), seq2.begin(), seq2.end());
  pack_set(both, codes, start, end, &num_seqs, &alen);
}

PackedAaSequenceSet::PackedAaSequenceSet(const std::vector<std::string> &seq)
  : PackedAaSequenceSet(sequence_views_from_strings(seq)) {}

PackedAaSequenceSet::PackedAaSequenceSet(
  const std::vector<std::string> &seq1,
  const std::vector<std::string> &seq2
) : PackedAaSequenceSet(
      sequence_views_from_strings(seq1),
      sequence_views_from_strings(seq2)
    ) {}

void PackedAaSequenceSet::fill_pos_major(
  int begin,
  int end,
  std::vector<uint8_t> &out
) const {
  const int n_j = end - begin;
  if (n_j <= 0 || alen <= 0) {
    out.clear();
    return;
  }
  out.assign(static_cast<std::size_t>(alen) * static_cast<std::size_t>(n_j), 0);
  for (int j = 0; j < n_j; ++j) {
    const uint8_t *src = row(begin + j);
    for (int k = 0; k < alen; ++k) {
      out[static_cast<std::size_t>(k) * static_cast<std::size_t>(n_j) +
          static_cast<std::size_t>(j)] = src[k];
    }
  }
}

void PackedAaSequenceSet::fill_pos_major_idx(
  const int *idx,
  int n_idx,
  std::vector<uint8_t> &out
) const {
  if (n_idx <= 0 || alen <= 0) {
    out.clear();
    return;
  }
  out.assign(static_cast<std::size_t>(alen) * static_cast<std::size_t>(n_idx), 0);
  for (int j = 0; j < n_idx; ++j) {
    const uint8_t *src = row(idx[j]);
    for (int k = 0; k < alen; ++k) {
      out[static_cast<std::size_t>(k) * static_cast<std::size_t>(n_idx) +
          static_cast<std::size_t>(j)] = src[k];
    }
  }
}

std::tuple<bool, int, double> PackedAaSequenceSet::success_score_and_dist(
  int i,
  int j,
  int min_overlap,
  bool ignore_gaps,
  const float *weight_lut
) const {
  const int lo = std::max(start[static_cast<std::size_t>(i)], start[static_cast<std::size_t>(j)]);
  const int hi = std::min(end[static_cast<std::size_t>(i)], end[static_cast<std::size_t>(j)]);
  const uint8_t *a = row(i);
  const uint8_t *b = row(j);
  int n_ok = 0;
  int n_match = 0;
  double acc = 0.0;
  int k = lo;
#if defined(__SSE2__)
  if (weight_lut == nullptr) {
    for (; k + 16 <= hi; k += 16) {
      identity_pair_16(a + k, b + k, ignore_gaps, &n_ok, &n_match);
    }
  }
#endif
  for (; k < hi; ++k) {
    if (weight_lut == nullptr) {
      identity_site(a[k], b[k], ignore_gaps, &n_ok, &n_match);
    } else {
      weighted_site(a[k], b[k], ignore_gaps, weight_lut, &n_ok, &acc);
    }
  }
  if (n_ok < min_overlap || n_ok == 0) {
    return {false, 0, 1.0};
  }
  const double d = weight_lut == nullptr
    ? 1.0 - static_cast<double>(n_match) / static_cast<double>(n_ok)
    : acc / static_cast<double>(n_ok);
  return {true, n_ok, d};
}

double PackedAaSequenceSet::dist(
  int i,
  int j,
  int min_overlap,
  bool ignore_gaps,
  const float *weight_lut
) const {
  return std::get<2>(
    success_score_and_dist(i, j, min_overlap, ignore_gaps, weight_lut)
  );
}

void PackedAaSequenceSet::score_query_vs_pos_major(
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
) const {
  if (n_target <= 0) {
    return;
  }
  const uint8_t *q = row(query_i);
  const int q_start = start[static_cast<std::size_t>(query_i)];
  const int q_end = end[static_cast<std::size_t>(query_i)];
  bool covered = ignore_gaps;
  if (!covered) {
    covered = true;
    for (int j = 0; j < n_target; ++j) {
      const int tj = target_idx[j];
      if (start[static_cast<std::size_t>(tj)] > q_start ||
          end[static_cast<std::size_t>(tj)] < q_end) {
        covered = false;
        break;
      }
    }
  }
  const bool use_vector = covered && alen <= 65535;

#if defined(__SSE2__)
  if (use_vector && weight_lut == nullptr) {
    thread_local std::vector<uint16_t> n_ok16;
    thread_local std::vector<uint16_t> n_match16;
    n_ok16.assign(static_cast<std::size_t>(n_target), 0);
    n_match16.assign(static_cast<std::size_t>(n_target), 0);
    for (int k = q_start; k < q_end; ++k) {
      const uint8_t a = q[k];
      if (ignore_gaps && a == 0) {
        continue;
      }
      const uint8_t *col = pos_major +
        static_cast<std::size_t>(k) * static_cast<std::size_t>(n_cols);
      int j = 0;
      for (; j + 16 <= n_target; j += 16) {
        identity_column_16(
          a,
          col + j,
          ignore_gaps,
          n_ok16.data() + j,
          n_match16.data() + j
        );
      }
      for (; j < n_target; ++j) {
        identity_site_u16(
          a,
          col[j],
          ignore_gaps,
          &n_ok16[static_cast<std::size_t>(j)],
          &n_match16[static_cast<std::size_t>(j)]
        );
      }
    }
    for (int j = 0; j < n_target; ++j) {
      const int nok = n_ok16[static_cast<std::size_t>(j)];
      n_ok_out[j] = nok;
      if (nok < min_overlap || nok == 0) {
        dist_out[j] = 1.0;
      } else {
        dist_out[j] = 1.0 -
          static_cast<double>(n_match16[static_cast<std::size_t>(j)]) /
          static_cast<double>(nok);
      }
    }
    return;
  }
#endif

#if defined(__AVX2__)
  if (use_vector && weight_lut != nullptr) {
    thread_local std::vector<int> n_ok;
    thread_local std::vector<double> acc;
    n_ok.assign(static_cast<std::size_t>(n_target), 0);
    acc.assign(static_cast<std::size_t>(n_target), 0.0);
    for (int k = q_start; k < q_end; ++k) {
      const uint8_t a = q[k];
      if (ignore_gaps && a == 0) {
        continue;
      }
      const uint8_t *col = pos_major +
        static_cast<std::size_t>(k) * static_cast<std::size_t>(n_cols);
      const float *prof_row = weight_lut + static_cast<std::size_t>(a) * kLut;
      int j = 0;
      for (; j + 8 <= n_target; j += 8) {
        weighted_column_8(
          prof_row,
          a,
          col + j,
          ignore_gaps,
          n_ok.data() + j,
          acc.data() + j
        );
      }
      for (; j < n_target; ++j) {
        const uint8_t b = col[j];
        if (a == 0 && b == 0) {
          continue;
        }
        if (ignore_gaps && (a == 0 || b == 0)) {
          continue;
        }
        ++n_ok[static_cast<std::size_t>(j)];
        acc[static_cast<std::size_t>(j)] += static_cast<double>(prof_row[b]);
      }
    }
    for (int j = 0; j < n_target; ++j) {
      const int nok = n_ok[static_cast<std::size_t>(j)];
      n_ok_out[j] = nok;
      if (nok < min_overlap || nok == 0) {
        dist_out[j] = 1.0;
      } else {
        dist_out[j] = acc[static_cast<std::size_t>(j)] /
          static_cast<double>(nok);
      }
    }
    return;
  }
#else
  (void)use_vector;
#endif

  std::vector<int> n_ok(static_cast<std::size_t>(n_target), 0);
  std::vector<int> n_match(static_cast<std::size_t>(n_target), 0);
  std::vector<double> acc(static_cast<std::size_t>(n_target), 0.0);
  std::vector<char> live(static_cast<std::size_t>(n_target), 1);
  int n_live = n_target;

  for (int k = 0; k < alen && n_live > 0; ++k) {
    if (k < q_start || k >= q_end) {
      continue;
    }
    const uint8_t a = q[k];
    const uint8_t *col =
      pos_major + static_cast<std::size_t>(k) * static_cast<std::size_t>(n_cols);
    const float *prof_row = nullptr;
    if (weight_lut != nullptr) {
      prof_row = weight_lut + static_cast<std::size_t>(a) * kLut;
    }
    for (int j = 0; j < n_target; ++j) {
      if (!live[static_cast<std::size_t>(j)]) {
        continue;
      }
      const int tj = target_idx[j];
      if (k < start[static_cast<std::size_t>(tj)] ||
          k >= end[static_cast<std::size_t>(tj)]) {
        continue;
      }
      const uint8_t b = col[j];
      if (weight_lut == nullptr) {
        identity_site(a, b, ignore_gaps, &n_ok[static_cast<std::size_t>(j)],
                      &n_match[static_cast<std::size_t>(j)]);
      } else {
        if (a == 0 && b == 0) {
          continue;
        }
        if (ignore_gaps && (a == 0 || b == 0)) {
          continue;
        }
        ++n_ok[static_cast<std::size_t>(j)];
        acc[static_cast<std::size_t>(j)] += static_cast<double>(prof_row[b]);
      }
    }
    if (early_exit && (k & 15) == 15) {
      const int remain = std::max(0, q_end - k - 1);
      for (int j = 0; j < n_target; ++j) {
        if (!live[static_cast<std::size_t>(j)]) {
          continue;
        }
        const int nok = n_ok[static_cast<std::size_t>(j)];
        if (nok == 0) {
          continue;
        }
        double best_d;
        if (weight_lut == nullptr) {
          best_d = 1.0 - static_cast<double>(n_match[static_cast<std::size_t>(j)] + remain) /
            static_cast<double>(nok + remain);
        } else {
          best_d = acc[static_cast<std::size_t>(j)] /
            static_cast<double>(nok + remain);
        }
        if (best_d > dist_threshold) {
          live[static_cast<std::size_t>(j)] = 0;
          --n_live;
        }
      }
    }
  }

  for (int j = 0; j < n_target; ++j) {
    const int nok = n_ok[static_cast<std::size_t>(j)];
    n_ok_out[j] = nok;
    if (nok < min_overlap || nok == 0) {
      dist_out[j] = 1.0;
    } else if (weight_lut == nullptr) {
      dist_out[j] = 1.0 - static_cast<double>(n_match[static_cast<std::size_t>(j)]) /
        static_cast<double>(nok);
    } else {
      dist_out[j] = acc[static_cast<std::size_t>(j)] / static_cast<double>(nok);
    }
  }
}
