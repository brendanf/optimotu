// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "hamming_stripe.h"
#include <algorithm>
#include <vector>

#if OPTIMOTU_HAMMING_KERNEL == 3 || OPTIMOTU_HAMMING_KERNEL == 2
#include <immintrin.h>
#endif
#if OPTIMOTU_HAMMING_KERNEL == 4
#include <arm_neon.h>
#endif

namespace {

void finalize_lane(
  uint32_t n_ok,
  uint32_t n_match,
  int min_overlap,
  double *dist_out,
  int *n_ok_out,
  int *n_match_out
) {
  *n_ok_out = static_cast<int>(n_ok);
  if (n_match_out != nullptr) {
    *n_match_out = static_cast<int>(n_match);
  }
  if (n_ok == 0 || static_cast<int>(n_ok) < min_overlap) {
    *dist_out = 1.0;
  } else {
    *dist_out = 1.0 - static_cast<double>(n_match) / static_cast<double>(n_ok);
  }
}

void score_stripe_scalar(
  const uint64_t *query_seq,
  const uint64_t *query_mask,
  const uint64_t *seq_major,
  const uint64_t *mask_major,
  int ulen,
  int mulen,
  int n_target,
  int n_cols,
  uint32_t *n_ok,
  uint32_t *n_match,
  int j0,
  int j1
) {
  for (int w = 0; w < mulen; ++w) {
    const uint64_t qm = query_mask[w];
    const uint64_t *row = mask_major +
      static_cast<std::size_t>(w) * static_cast<std::size_t>(n_cols);
    for (int j = j0; j < j1; ++j) {
      n_ok[j] += static_cast<uint32_t>(__builtin_popcountll(qm & row[j]));
    }
  }
  for (int w = 0; w < ulen; ++w) {
    const uint64_t qs = query_seq[w];
    const uint64_t *row = seq_major +
      static_cast<std::size_t>(w) * static_cast<std::size_t>(n_cols);
    for (int j = j0; j < j1; ++j) {
      n_match[j] += static_cast<uint32_t>(__builtin_popcountll(qs & row[j]));
    }
  }
}

#if OPTIMOTU_HAMMING_KERNEL == 2
void score_stripe_avx2(
  const uint64_t *query_seq,
  const uint64_t *query_mask,
  const uint64_t *seq_major,
  const uint64_t *mask_major,
  int ulen,
  int mulen,
  int n_target,
  int n_cols,
  uint32_t *n_ok,
  uint32_t *n_match
) {
  int j = 0;
  for (; j + 4 <= n_target; j += 4) {
    for (int w = 0; w < mulen; ++w) {
      const __m256i qm = _mm256_set1_epi64x(
        static_cast<long long>(query_mask[w])
      );
      const uint64_t *row = mask_major +
        static_cast<std::size_t>(w) * static_cast<std::size_t>(n_cols) +
        static_cast<std::size_t>(j);
      const __m256i tm = _mm256_loadu_si256(
        reinterpret_cast<const __m256i *>(row)
      );
      const __m256i vand = _mm256_and_si256(qm, tm);
      n_ok[j] += static_cast<uint32_t>(
        __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 0)))
      );
      n_ok[j + 1] += static_cast<uint32_t>(
        __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 1)))
      );
      n_ok[j + 2] += static_cast<uint32_t>(
        __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 2)))
      );
      n_ok[j + 3] += static_cast<uint32_t>(
        __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 3)))
      );
    }
    for (int w = 0; w < ulen; ++w) {
      const __m256i qs = _mm256_set1_epi64x(
        static_cast<long long>(query_seq[w])
      );
      const uint64_t *row = seq_major +
        static_cast<std::size_t>(w) * static_cast<std::size_t>(n_cols) +
        static_cast<std::size_t>(j);
      const __m256i ts = _mm256_loadu_si256(
        reinterpret_cast<const __m256i *>(row)
      );
      const __m256i vand = _mm256_and_si256(qs, ts);
      n_match[j] += static_cast<uint32_t>(
        __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 0)))
      );
      n_match[j + 1] += static_cast<uint32_t>(
        __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 1)))
      );
      n_match[j + 2] += static_cast<uint32_t>(
        __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 2)))
      );
      n_match[j + 3] += static_cast<uint32_t>(
        __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 3)))
      );
    }
  }
  if (j < n_target) {
    score_stripe_scalar(
      query_seq, query_mask, seq_major, mask_major, ulen, mulen,
      n_target, n_cols, n_ok, n_match, j, n_target
    );
  }
}
#endif

#if OPTIMOTU_HAMMING_KERNEL == 3
void score_stripe_avx512(
  const uint64_t *query_seq,
  const uint64_t *query_mask,
  const uint64_t *seq_major,
  const uint64_t *mask_major,
  int ulen,
  int mulen,
  int n_target,
  int n_cols,
  uint32_t *n_ok,
  uint32_t *n_match
) {
  int j = 0;
  for (; j + 8 <= n_target; j += 8) {
    __m256i ok = _mm256_loadu_si256(
      reinterpret_cast<const __m256i *>(n_ok + j)
    );
    for (int w = 0; w < mulen; ++w) {
      const __m512i qm = _mm512_set1_epi64(
        static_cast<long long>(query_mask[w])
      );
      const uint64_t *row = mask_major +
        static_cast<std::size_t>(w) * static_cast<std::size_t>(n_cols) +
        static_cast<std::size_t>(j);
      const __m512i tm = _mm512_loadu_si512(row);
      const __m512i cnt = _mm512_popcnt_epi64(_mm512_and_si512(qm, tm));
      ok = _mm256_add_epi32(ok, _mm512_cvtepi64_epi32(cnt));
    }
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(n_ok + j), ok);

    __m256i mt = _mm256_loadu_si256(
      reinterpret_cast<const __m256i *>(n_match + j)
    );
    for (int w = 0; w < ulen; ++w) {
      const __m512i qs = _mm512_set1_epi64(
        static_cast<long long>(query_seq[w])
      );
      const uint64_t *row = seq_major +
        static_cast<std::size_t>(w) * static_cast<std::size_t>(n_cols) +
        static_cast<std::size_t>(j);
      const __m512i ts = _mm512_loadu_si512(row);
      const __m512i cnt = _mm512_popcnt_epi64(_mm512_and_si512(qs, ts));
      mt = _mm256_add_epi32(mt, _mm512_cvtepi64_epi32(cnt));
    }
    _mm256_storeu_si256(reinterpret_cast<__m256i *>(n_match + j), mt);
  }
  if (j < n_target) {
    score_stripe_scalar(
      query_seq, query_mask, seq_major, mask_major, ulen, mulen,
      n_target, n_cols, n_ok, n_match, j, n_target
    );
  }
}
#endif

#if OPTIMOTU_HAMMING_KERNEL == 4
void neon_popcnt2(uint64x2_t v, uint32_t *a, uint32_t *b) {
  const uint8x16_t c = vcntq_u8(vreinterpretq_u8_u64(v));
  const uint16x8_t s1 = vpaddlq_u8(c);
  const uint32x4_t s2 = vpaddlq_u16(s1);
  const uint64x2_t s3 = vpaddlq_u32(s2);
  *a += static_cast<uint32_t>(vgetq_lane_u64(s3, 0));
  *b += static_cast<uint32_t>(vgetq_lane_u64(s3, 1));
}

void score_stripe_neon(
  const uint64_t *query_seq,
  const uint64_t *query_mask,
  const uint64_t *seq_major,
  const uint64_t *mask_major,
  int ulen,
  int mulen,
  int n_target,
  int n_cols,
  uint32_t *n_ok,
  uint32_t *n_match
) {
  int j = 0;
  for (; j + 2 <= n_target; j += 2) {
    for (int w = 0; w < mulen; ++w) {
      const uint64x2_t qm = vdupq_n_u64(query_mask[w]);
      const uint64_t *row = mask_major +
        static_cast<std::size_t>(w) * static_cast<std::size_t>(n_cols) +
        static_cast<std::size_t>(j);
      const uint64x2_t tm = vld1q_u64(row);
      neon_popcnt2(vandq_u64(qm, tm), &n_ok[j], &n_ok[j + 1]);
    }
    for (int w = 0; w < ulen; ++w) {
      const uint64x2_t qs = vdupq_n_u64(query_seq[w]);
      const uint64_t *row = seq_major +
        static_cast<std::size_t>(w) * static_cast<std::size_t>(n_cols) +
        static_cast<std::size_t>(j);
      const uint64x2_t ts = vld1q_u64(row);
      neon_popcnt2(vandq_u64(qs, ts), &n_match[j], &n_match[j + 1]);
    }
  }
  if (j < n_target) {
    score_stripe_scalar(
      query_seq, query_mask, seq_major, mask_major, ulen, mulen,
      n_target, n_cols, n_ok, n_match, j, n_target
    );
  }
}
#endif

}  // namespace

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
) {
  if (n_target <= 0) {
    return;
  }
  std::vector<uint32_t> n_ok(static_cast<std::size_t>(n_target), 0);
  std::vector<uint32_t> n_match(static_cast<std::size_t>(n_target), 0);
#if OPTIMOTU_HAMMING_KERNEL == 3
  score_stripe_avx512(
    query_seq, query_mask, seq_major, mask_major, ulen, mulen,
    n_target, n_cols, n_ok.data(), n_match.data()
  );
#elif OPTIMOTU_HAMMING_KERNEL == 2
  score_stripe_avx2(
    query_seq, query_mask, seq_major, mask_major, ulen, mulen,
    n_target, n_cols, n_ok.data(), n_match.data()
  );
#elif OPTIMOTU_HAMMING_KERNEL == 4
  score_stripe_neon(
    query_seq, query_mask, seq_major, mask_major, ulen, mulen,
    n_target, n_cols, n_ok.data(), n_match.data()
  );
#else
  score_stripe_scalar(
    query_seq, query_mask, seq_major, mask_major, ulen, mulen,
    n_target, n_cols, n_ok.data(), n_match.data(), 0, n_target
  );
#endif
  for (int j = 0; j < n_target; ++j) {
    finalize_lane(
      n_ok[static_cast<std::size_t>(j)],
      n_match[static_cast<std::size_t>(j)],
      min_overlap,
      dist_out + j,
      n_ok_out + j,
      n_match_out == nullptr ? nullptr : n_match_out + j
    );
  }
}

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
) {
  uint32_t n_ok = 0;
  uint32_t n_match = 0;
  int w = 0;
#if OPTIMOTU_HAMMING_KERNEL == 3
  for (; w + 8 <= mulen; w += 8) {
    const __m512i va = _mm512_loadu_si512(ma + w);
    const __m512i vb = _mm512_loadu_si512(mb + w);
    const __m512i cnt = _mm512_popcnt_epi64(_mm512_and_si512(va, vb));
    alignas(64) uint64_t tmp[8];
    _mm512_store_si512(tmp, cnt);
    for (int k = 0; k < 8; ++k) {
      n_ok += static_cast<uint32_t>(tmp[k]);
    }
  }
#elif OPTIMOTU_HAMMING_KERNEL == 2
  for (; w + 4 <= mulen; w += 4) {
    const __m256i va = _mm256_loadu_si256(
      reinterpret_cast<const __m256i *>(ma + w)
    );
    const __m256i vb = _mm256_loadu_si256(
      reinterpret_cast<const __m256i *>(mb + w)
    );
    const __m256i vand = _mm256_and_si256(va, vb);
    n_ok += static_cast<uint32_t>(
      __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 0)))
    );
    n_ok += static_cast<uint32_t>(
      __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 1)))
    );
    n_ok += static_cast<uint32_t>(
      __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 2)))
    );
    n_ok += static_cast<uint32_t>(
      __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 3)))
    );
  }
#elif OPTIMOTU_HAMMING_KERNEL == 4
  for (; w + 2 <= mulen; w += 2) {
    const uint64x2_t va = vld1q_u64(ma + w);
    const uint64x2_t vb = vld1q_u64(mb + w);
    uint32_t x = 0;
    uint32_t y = 0;
    neon_popcnt2(vandq_u64(va, vb), &x, &y);
    n_ok += x + y;
  }
#endif
  for (; w < mulen; ++w) {
    n_ok += static_cast<uint32_t>(__builtin_popcountll(ma[w] & mb[w]));
  }

  w = 0;
#if OPTIMOTU_HAMMING_KERNEL == 3
  for (; w + 8 <= ulen; w += 8) {
    const __m512i va = _mm512_loadu_si512(a + w);
    const __m512i vb = _mm512_loadu_si512(b + w);
    const __m512i cnt = _mm512_popcnt_epi64(_mm512_and_si512(va, vb));
    alignas(64) uint64_t tmp[8];
    _mm512_store_si512(tmp, cnt);
    for (int k = 0; k < 8; ++k) {
      n_match += static_cast<uint32_t>(tmp[k]);
    }
  }
#elif OPTIMOTU_HAMMING_KERNEL == 2
  for (; w + 4 <= ulen; w += 4) {
    const __m256i va = _mm256_loadu_si256(
      reinterpret_cast<const __m256i *>(a + w)
    );
    const __m256i vb = _mm256_loadu_si256(
      reinterpret_cast<const __m256i *>(b + w)
    );
    const __m256i vand = _mm256_and_si256(va, vb);
    n_match += static_cast<uint32_t>(
      __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 0)))
    );
    n_match += static_cast<uint32_t>(
      __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 1)))
    );
    n_match += static_cast<uint32_t>(
      __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 2)))
    );
    n_match += static_cast<uint32_t>(
      __builtin_popcountll(static_cast<uint64_t>(_mm256_extract_epi64(vand, 3)))
    );
  }
#elif OPTIMOTU_HAMMING_KERNEL == 4
  for (; w + 2 <= ulen; w += 2) {
    const uint64x2_t va = vld1q_u64(a + w);
    const uint64x2_t vb = vld1q_u64(b + w);
    uint32_t x = 0;
    uint32_t y = 0;
    neon_popcnt2(vandq_u64(va, vb), &x, &y);
    n_match += x + y;
  }
#endif
  for (; w < ulen; ++w) {
    n_match += static_cast<uint32_t>(__builtin_popcountll(a[w] & b[w]));
  }

  *num_ok = static_cast<int>(n_ok);
  *num_matches = static_cast<int>(n_match);
  if (n_ok == 0 || static_cast<int>(n_ok) < min_overlap) {
    *dist = 1.0;
    return 1;
  }
  *dist = 1.0 - static_cast<double>(n_match) / static_cast<double>(n_ok);
  return 0;
}
