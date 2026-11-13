// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "PackedAaSequenceSet.h"
#include "aa_kernel.h"
#include "optimotu.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

#if OPTIMOTU_AA_KERNEL >= 3
#include <immintrin.h>
#elif OPTIMOTU_AA_KERNEL == 2
#if defined(__ARM_NEON)
#include <arm_neon.h>
#else
#include <emmintrin.h>
#endif
#endif

namespace {

inline bool evo_site_ok(const EvoModel &model, uint8_t a, uint8_t b) {
  return a < EvoModel::kCode && b < EvoModel::kCode &&
    model.is_state[a] != 0 && model.is_state[b] != 0;
}

inline int overlap_lo(const PackedAaSequenceSet &pss, int i, int j) {
  return std::max(
    pss.start[static_cast<std::size_t>(i)],
    pss.start[static_cast<std::size_t>(j)]
  );
}

inline int overlap_hi(const PackedAaSequenceSet &pss, int i, int j) {
  return std::min(
    pss.end[static_cast<std::size_t>(i)],
    pss.end[static_cast<std::size_t>(j)]
  );
}

void fill_counts(
  const PackedAaSequenceSet &pss,
  int i,
  int j,
  const EvoModel &model,
  int counts[EvoModel::kState][EvoModel::kState]
) {
  const uint8_t *a = pss.row(i);
  const uint8_t *b = pss.row(j);
  const int lo = overlap_lo(pss, i, j);
  const int hi = overlap_hi(pss, i, j);
  for (int k = lo; k < hi; ++k) {
    const uint8_t ca = a[k];
    const uint8_t cb = b[k];
    if (!evo_site_ok(model, ca, cb)) {
      continue;
    }
    const int si = model.code_to_state[ca];
    const int sj = model.code_to_state[cb];
    ++counts[si][sj];
  }
}

// Batch width for exp(lambda * t). Dots over the 20 eigenvalues use the
// same ISA, but they run inside one pair.
#if OPTIMOTU_AA_KERNEL == 4
constexpr int kNewtonBatch = 8;
#elif OPTIMOTU_AA_KERNEL == 3
constexpr int kNewtonBatch = 4;
#elif OPTIMOTU_AA_KERNEL == 2
constexpr int kNewtonBatch = 2;
#else
constexpr int kNewtonBatch = 1;
#endif

constexpr double kLn2Hi = 0.69314718055994528622676398299518041312694549560546875;
constexpr double kLn2Lo = 2.319046813846299558417771068559111445533847808837890625e-17;
constexpr double kInvLn2 = 1.44269504088896340735992468100189213771759033203125;

// V[i][k] and U[j][k] are contiguous in the eigenvalue index.
struct SpectralRows {
  double v[EvoModel::kState][EvoModel::kState];
  double u[EvoModel::kState][EvoModel::kState];
};

struct NewtonHit {
  int i;
  int j;
  int c;
};

struct NewtonLane {
  double t = 0.0;
  double t_lo = 0.0;
  double t_hi = 0.0;
  double grid_t = 0.0;
  double grid_ll = 0.0;
  double out = 0.0;
  bool present = false;
  bool moving = false;
  int n_hit = 0;
  NewtonHit hits[EvoModel::kState * EvoModel::kState];
  int counts[EvoModel::kState][EvoModel::kState] = {};
};

SpectralRows spectral_rows(const EvoModel &model) {
  SpectralRows rows;
  for (int i = 0; i < EvoModel::kState; ++i) {
    for (int k = 0; k < EvoModel::kState; ++k) {
      rows.v[i][k] = model.V_ik(i, k);
      rows.u[i][k] = model.Vinv_kj(k, i);
    }
  }
  return rows;
}

void collect_hits(NewtonLane &lane) {
  lane.n_hit = 0;
  for (int i = 0; i < EvoModel::kState; ++i) {
    for (int j = 0; j < EvoModel::kState; ++j) {
      const int c = lane.counts[i][j];
      if (c == 0) {
        continue;
      }
      lane.hits[lane.n_hit++] = NewtonHit{i, j, c};
    }
  }
}

inline double dot_eigen(
  const double *a,
  const double *b,
  const double *c
) {
#if OPTIMOTU_AA_KERNEL == 4
  __m512d acc = _mm512_setzero_pd();
  for (int k = 0; k < 16; k += 8) {
    const __m512d va = _mm512_loadu_pd(a + k);
    const __m512d vb = _mm512_loadu_pd(b + k);
    const __m512d vc = _mm512_loadu_pd(c + k);
    acc = _mm512_add_pd(acc, _mm512_mul_pd(va, _mm512_mul_pd(vb, vc)));
  }
  __m256d tail = _mm256_mul_pd(
    _mm256_loadu_pd(a + 16),
    _mm256_mul_pd(_mm256_loadu_pd(b + 16), _mm256_loadu_pd(c + 16))
  );
  const __m256d hi = _mm512_extractf64x4_pd(acc, 1);
  __m256d sum = _mm256_add_pd(_mm512_castpd512_pd256(acc), hi);
  sum = _mm256_add_pd(sum, tail);
  const __m128d lo = _mm_add_pd(
    _mm256_castpd256_pd128(sum),
    _mm256_extractf128_pd(sum, 1)
  );
  const __m128d folded = _mm_add_sd(lo, _mm_unpackhi_pd(lo, lo));
  return _mm_cvtsd_f64(folded);
#elif OPTIMOTU_AA_KERNEL == 3
  __m256d acc = _mm256_setzero_pd();
  for (int k = 0; k < EvoModel::kState; k += 4) {
    const __m256d va = _mm256_loadu_pd(a + k);
    const __m256d vb = _mm256_loadu_pd(b + k);
    const __m256d vc = _mm256_loadu_pd(c + k);
    acc = _mm256_add_pd(acc, _mm256_mul_pd(va, _mm256_mul_pd(vb, vc)));
  }
  const __m128d lo = _mm_add_pd(
    _mm256_castpd256_pd128(acc),
    _mm256_extractf128_pd(acc, 1)
  );
  const __m128d folded = _mm_add_sd(lo, _mm_unpackhi_pd(lo, lo));
  return _mm_cvtsd_f64(folded);
#elif OPTIMOTU_AA_KERNEL == 2 && defined(__ARM_NEON)
  float64x2_t acc = vdupq_n_f64(0.0);
  for (int k = 0; k < EvoModel::kState; k += 2) {
    const float64x2_t va = vld1q_f64(a + k);
    const float64x2_t vb = vld1q_f64(b + k);
    const float64x2_t vc = vld1q_f64(c + k);
    acc = vaddq_f64(acc, vmulq_f64(va, vmulq_f64(vb, vc)));
  }
  return vgetq_lane_f64(acc, 0) + vgetq_lane_f64(acc, 1);
#elif OPTIMOTU_AA_KERNEL == 2
  __m128d acc = _mm_setzero_pd();
  for (int k = 0; k < EvoModel::kState; k += 2) {
    const __m128d va = _mm_loadu_pd(a + k);
    const __m128d vb = _mm_loadu_pd(b + k);
    const __m128d vc = _mm_loadu_pd(c + k);
    acc = _mm_add_pd(acc, _mm_mul_pd(va, _mm_mul_pd(vb, vc)));
  }
  const __m128d folded = _mm_add_sd(acc, _mm_unpackhi_pd(acc, acc));
  return _mm_cvtsd_f64(folded);
#else
  double sum = 0.0;
  for (int k = 0; k < EvoModel::kState; ++k) {
    sum += a[k] * b[k] * c[k];
  }
  return sum;
#endif
}

// p, dp, and d2p share V[i] * U[j], so one pass computes all three.
inline void dot_eigen3(
  const double *v,
  const double *u,
  const double *e,
  const double *de,
  const double *d2e,
  double *p_out,
  double *dp_out,
  double *d2p_out
) {
#if OPTIMOTU_AA_KERNEL == 4
  __m512d ap = _mm512_setzero_pd();
  __m512d adp = _mm512_setzero_pd();
  __m512d ad2 = _mm512_setzero_pd();
  for (int k = 0; k < 16; k += 8) {
    const __m512d w = _mm512_mul_pd(
      _mm512_loadu_pd(v + k),
      _mm512_loadu_pd(u + k)
    );
    ap = _mm512_add_pd(ap, _mm512_mul_pd(w, _mm512_loadu_pd(e + k)));
    adp = _mm512_add_pd(adp, _mm512_mul_pd(w, _mm512_loadu_pd(de + k)));
    ad2 = _mm512_add_pd(ad2, _mm512_mul_pd(w, _mm512_loadu_pd(d2e + k)));
  }
  const __m256d w4 = _mm256_mul_pd(
    _mm256_loadu_pd(v + 16),
    _mm256_loadu_pd(u + 16)
  );
  const __m256d tp = _mm256_mul_pd(w4, _mm256_loadu_pd(e + 16));
  const __m256d tdp = _mm256_mul_pd(w4, _mm256_loadu_pd(de + 16));
  const __m256d td2 = _mm256_mul_pd(w4, _mm256_loadu_pd(d2e + 16));
  auto fold = [](__m512d acc, __m256d tail) {
    const __m256d hi = _mm512_extractf64x4_pd(acc, 1);
    __m256d sum = _mm256_add_pd(_mm512_castpd512_pd256(acc), hi);
    sum = _mm256_add_pd(sum, tail);
    const __m128d lo = _mm_add_pd(
      _mm256_castpd256_pd128(sum),
      _mm256_extractf128_pd(sum, 1)
    );
    return _mm_cvtsd_f64(_mm_add_sd(lo, _mm_unpackhi_pd(lo, lo)));
  };
  *p_out = fold(ap, tp);
  *dp_out = fold(adp, tdp);
  *d2p_out = fold(ad2, td2);
#elif OPTIMOTU_AA_KERNEL == 3
  __m256d ap = _mm256_setzero_pd();
  __m256d adp = _mm256_setzero_pd();
  __m256d ad2 = _mm256_setzero_pd();
  for (int k = 0; k < EvoModel::kState; k += 4) {
    const __m256d w = _mm256_mul_pd(
      _mm256_loadu_pd(v + k),
      _mm256_loadu_pd(u + k)
    );
    ap = _mm256_add_pd(ap, _mm256_mul_pd(w, _mm256_loadu_pd(e + k)));
    adp = _mm256_add_pd(adp, _mm256_mul_pd(w, _mm256_loadu_pd(de + k)));
    ad2 = _mm256_add_pd(ad2, _mm256_mul_pd(w, _mm256_loadu_pd(d2e + k)));
  }
  auto fold = [](__m256d acc) {
    const __m128d lo = _mm_add_pd(
      _mm256_castpd256_pd128(acc),
      _mm256_extractf128_pd(acc, 1)
    );
    return _mm_cvtsd_f64(_mm_add_sd(lo, _mm_unpackhi_pd(lo, lo)));
  };
  *p_out = fold(ap);
  *dp_out = fold(adp);
  *d2p_out = fold(ad2);
#elif OPTIMOTU_AA_KERNEL == 2 && defined(__ARM_NEON)
  float64x2_t ap = vdupq_n_f64(0.0);
  float64x2_t adp = vdupq_n_f64(0.0);
  float64x2_t ad2 = vdupq_n_f64(0.0);
  for (int k = 0; k < EvoModel::kState; k += 2) {
    const float64x2_t w = vmulq_f64(vld1q_f64(v + k), vld1q_f64(u + k));
    ap = vaddq_f64(ap, vmulq_f64(w, vld1q_f64(e + k)));
    adp = vaddq_f64(adp, vmulq_f64(w, vld1q_f64(de + k)));
    ad2 = vaddq_f64(ad2, vmulq_f64(w, vld1q_f64(d2e + k)));
  }
  auto fold = [](float64x2_t acc) {
    return vgetq_lane_f64(acc, 0) + vgetq_lane_f64(acc, 1);
  };
  *p_out = fold(ap);
  *dp_out = fold(adp);
  *d2p_out = fold(ad2);
#elif OPTIMOTU_AA_KERNEL == 2
  __m128d ap = _mm_setzero_pd();
  __m128d adp = _mm_setzero_pd();
  __m128d ad2 = _mm_setzero_pd();
  for (int k = 0; k < EvoModel::kState; k += 2) {
    const __m128d w = _mm_mul_pd(_mm_loadu_pd(v + k), _mm_loadu_pd(u + k));
    ap = _mm_add_pd(ap, _mm_mul_pd(w, _mm_loadu_pd(e + k)));
    adp = _mm_add_pd(adp, _mm_mul_pd(w, _mm_loadu_pd(de + k)));
    ad2 = _mm_add_pd(ad2, _mm_mul_pd(w, _mm_loadu_pd(d2e + k)));
  }
  auto fold = [](__m128d acc) {
    return _mm_cvtsd_f64(_mm_add_sd(acc, _mm_unpackhi_pd(acc, acc)));
  };
  *p_out = fold(ap);
  *dp_out = fold(adp);
  *d2p_out = fold(ad2);
#else
  double p = 0.0;
  double dp = 0.0;
  double d2p = 0.0;
  for (int k = 0; k < EvoModel::kState; ++k) {
    const double w = v[k] * u[k];
    p += w * e[k];
    dp += w * de[k];
    d2p += w * d2e[k];
  }
  *p_out = p;
  *dp_out = dp;
  *d2p_out = d2p;
#endif
}

#if OPTIMOTU_AA_KERNEL == 4

inline __m512d horner_exp_8(__m512d f) {
  __m512d p = _mm512_set1_pd(1.0 / 479001600.0);
  const double c[] = {
    1.0 / 39916800.0,
    1.0 / 3628800.0,
    1.0 / 362880.0,
    1.0 / 40320.0,
    1.0 / 5040.0,
    1.0 / 720.0,
    1.0 / 120.0,
    1.0 / 24.0,
    1.0 / 6.0,
    0.5,
    1.0,
    1.0
  };
  for (double coeff : c) {
    p = _mm512_add_pd(_mm512_mul_pd(p, f), _mm512_set1_pd(coeff));
  }
  return p;
}

inline __m512d exp_batch_8(__m512d x) {
  const __mmask8 too_hi = _mm512_cmp_pd_mask(x, _mm512_set1_pd(709.0), _CMP_GT_OQ);
  const __mmask8 too_lo = _mm512_cmp_pd_mask(x, _mm512_set1_pd(-745.0), _CMP_LT_OQ);
  const __m512d n = _mm512_roundscale_pd(
    _mm512_mul_pd(x, _mm512_set1_pd(kInvLn2)),
    _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC
  );
  __m512d f = _mm512_sub_pd(x, _mm512_mul_pd(n, _mm512_set1_pd(kLn2Hi)));
  f = _mm512_sub_pd(f, _mm512_mul_pd(n, _mm512_set1_pd(kLn2Lo)));
  const __m512d p = horner_exp_8(f);
  __m512i bits = _mm512_cvtepi32_epi64(_mm512_cvtpd_epi32(n));
  bits = _mm512_add_epi64(bits, _mm512_set1_epi64(1023));
  bits = _mm512_slli_epi64(bits, 52);
  __m512d y = _mm512_mul_pd(p, _mm512_castsi512_pd(bits));
  y = _mm512_mask_blend_pd(too_hi, y, _mm512_set1_pd(std::numeric_limits<double>::infinity()));
  return _mm512_mask_blend_pd(too_lo, y, _mm512_setzero_pd());
}

#elif OPTIMOTU_AA_KERNEL == 3

inline __m256d horner_exp_4(__m256d f) {
  __m256d p = _mm256_set1_pd(1.0 / 479001600.0);
  const double c[] = {
    1.0 / 39916800.0,
    1.0 / 3628800.0,
    1.0 / 362880.0,
    1.0 / 40320.0,
    1.0 / 5040.0,
    1.0 / 720.0,
    1.0 / 120.0,
    1.0 / 24.0,
    1.0 / 6.0,
    0.5,
    1.0,
    1.0
  };
  for (double coeff : c) {
    p = _mm256_add_pd(_mm256_mul_pd(p, f), _mm256_set1_pd(coeff));
  }
  return p;
}

inline __m256d exp_batch_4(__m256d x) {
  const __m256d too_hi = _mm256_cmp_pd(x, _mm256_set1_pd(709.0), _CMP_GT_OQ);
  const __m256d too_lo = _mm256_cmp_pd(x, _mm256_set1_pd(-745.0), _CMP_LT_OQ);
  const __m256d n = _mm256_round_pd(
    _mm256_mul_pd(x, _mm256_set1_pd(kInvLn2)),
    _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC
  );
  __m256d f = _mm256_sub_pd(x, _mm256_mul_pd(n, _mm256_set1_pd(kLn2Hi)));
  f = _mm256_sub_pd(f, _mm256_mul_pd(n, _mm256_set1_pd(kLn2Lo)));
  const __m256d p = horner_exp_4(f);
  const __m128i ni = _mm256_cvtpd_epi32(n);
  __m256i bits = _mm256_cvtepi32_epi64(ni);
  bits = _mm256_add_epi64(bits, _mm256_set1_epi64x(1023));
  bits = _mm256_slli_epi64(bits, 52);
  __m256d y = _mm256_mul_pd(p, _mm256_castsi256_pd(bits));
  const __m256d inf = _mm256_set1_pd(std::numeric_limits<double>::infinity());
  y = _mm256_blendv_pd(y, inf, too_hi);
  return _mm256_blendv_pd(y, _mm256_setzero_pd(), too_lo);
}

#elif OPTIMOTU_AA_KERNEL == 2 && defined(__ARM_NEON)

inline float64x2_t horner_exp_2(float64x2_t f) {
  float64x2_t p = vdupq_n_f64(1.0 / 479001600.0);
  const double c[] = {
    1.0 / 39916800.0,
    1.0 / 3628800.0,
    1.0 / 362880.0,
    1.0 / 40320.0,
    1.0 / 5040.0,
    1.0 / 720.0,
    1.0 / 120.0,
    1.0 / 24.0,
    1.0 / 6.0,
    0.5,
    1.0,
    1.0
  };
  for (double coeff : c) {
    p = vaddq_f64(vmulq_f64(p, f), vdupq_n_f64(coeff));
  }
  return p;
}

inline float64x2_t exp_batch_2(float64x2_t x) {
  const float64x2_t n = vrndnq_f64(vmulq_f64(x, vdupq_n_f64(kInvLn2)));
  float64x2_t f = vsubq_f64(x, vmulq_f64(n, vdupq_n_f64(kLn2Hi)));
  f = vsubq_f64(f, vmulq_f64(n, vdupq_n_f64(kLn2Lo)));
  const float64x2_t p = horner_exp_2(f);
  int64x2_t bits = vcvtq_s64_f64(n);
  bits = vaddq_s64(bits, vdupq_n_s64(1023));
  bits = vshlq_n_s64(bits, 52);
  return vmulq_f64(p, vreinterpretq_f64_s64(bits));
}

#elif OPTIMOTU_AA_KERNEL == 2

inline __m128d horner_exp_2(__m128d f) {
  __m128d p = _mm_set1_pd(1.0 / 479001600.0);
  const double c[] = {
    1.0 / 39916800.0,
    1.0 / 3628800.0,
    1.0 / 362880.0,
    1.0 / 40320.0,
    1.0 / 5040.0,
    1.0 / 720.0,
    1.0 / 120.0,
    1.0 / 24.0,
    1.0 / 6.0,
    0.5,
    1.0,
    1.0
  };
  for (double coeff : c) {
    p = _mm_add_pd(_mm_mul_pd(p, f), _mm_set1_pd(coeff));
  }
  return p;
}

inline __m128d sse2_round_nearest(__m128d x) {
  const __m128d sign = _mm_and_pd(x, _mm_set1_pd(-0.0));
  const __m128d magic = _mm_or_pd(_mm_set1_pd(4503599627370496.0), sign);
  return _mm_sub_pd(_mm_add_pd(x, magic), magic);
}

inline double scale_pow2(int exponent) {
  if (exponent > 1023) {
    return std::numeric_limits<double>::infinity();
  }
  if (exponent < -1022) {
    return 0.0;
  }
  const std::uint64_t bits =
    static_cast<std::uint64_t>(exponent + 1023) << 52;
  double out = 0.0;
  std::memcpy(&out, &bits, sizeof(out));
  return out;
}

inline __m128d exp_batch_2(__m128d x) {
  const __m128d n = sse2_round_nearest(_mm_mul_pd(x, _mm_set1_pd(kInvLn2)));
  __m128d f = _mm_sub_pd(x, _mm_mul_pd(n, _mm_set1_pd(kLn2Hi)));
  f = _mm_sub_pd(f, _mm_mul_pd(n, _mm_set1_pd(kLn2Lo)));
  const __m128d p = horner_exp_2(f);
  const __m128i ni = _mm_cvttpd_epi32(n);
  const int n0 = _mm_cvtsi128_si32(ni);
  const int n1 = _mm_cvtsi128_si32(_mm_srli_si128(ni, 4));
  const __m128d scale = _mm_set_pd(scale_pow2(n1), scale_pow2(n0));
  __m128d y = _mm_mul_pd(p, scale);
  const __m128d too_hi = _mm_cmpgt_pd(x, _mm_set1_pd(709.0));
  const __m128d too_lo = _mm_cmplt_pd(x, _mm_set1_pd(-745.0));
  const __m128d inf = _mm_set1_pd(std::numeric_limits<double>::infinity());
  y = _mm_or_pd(_mm_andnot_pd(too_hi, y), _mm_and_pd(too_hi, inf));
  return _mm_andnot_pd(too_lo, y);
}

#endif

// One SIMD exp per eigenvalue, across the batch. Each pair then reuses
// its own 20 values for every occupied residue pair.
void hoist_exponentials(
  const NewtonLane *lanes,
  const EvoModel &model,
  double e[kNewtonBatch][EvoModel::kState],
  double de[kNewtonBatch][EvoModel::kState],
  double d2e[kNewtonBatch][EvoModel::kState]
) {
  alignas(64) double tpack[kNewtonBatch];
  for (int lane = 0; lane < kNewtonBatch; ++lane) {
    tpack[lane] = lanes[lane].present ? lanes[lane].t : 0.0;
  }
  for (int k = 0; k < EvoModel::kState; ++k) {
    const double lam = model.lambda[static_cast<std::size_t>(k)];
    alignas(64) double ek[kNewtonBatch];
#if OPTIMOTU_AA_KERNEL == 4
    const __m512d y = exp_batch_8(
      _mm512_mul_pd(_mm512_load_pd(tpack), _mm512_set1_pd(lam))
    );
    _mm512_store_pd(ek, y);
#elif OPTIMOTU_AA_KERNEL == 3
    const __m256d y = exp_batch_4(
      _mm256_mul_pd(_mm256_load_pd(tpack), _mm256_set1_pd(lam))
    );
    _mm256_store_pd(ek, y);
#elif OPTIMOTU_AA_KERNEL == 2 && defined(__ARM_NEON)
    const float64x2_t y = exp_batch_2(
      vmulq_f64(vld1q_f64(tpack), vdupq_n_f64(lam))
    );
    vst1q_f64(ek, y);
#elif OPTIMOTU_AA_KERNEL == 2
    const __m128d y = exp_batch_2(
      _mm_mul_pd(_mm_load_pd(tpack), _mm_set1_pd(lam))
    );
    _mm_store_pd(ek, y);
#else
    ek[0] = std::exp(lam * tpack[0]);
#endif
    for (int lane = 0; lane < kNewtonBatch; ++lane) {
      e[lane][k] = ek[lane];
      de[lane][k] = ek[lane] * lam;
      d2e[lane][k] = de[lane][k] * lam;
    }
  }
}

double lane_loglik(
  const NewtonLane &lane,
  const SpectralRows &rows,
  const double *e
) {
  double sum = 0.0;
  for (int h = 0; h < lane.n_hit; ++h) {
    const NewtonHit &hit = lane.hits[h];
    double p = dot_eigen(rows.v[hit.i], rows.u[hit.j], e);
    if (!(p > 1e-300)) {
      p = 1e-300;
    }
    sum += static_cast<double>(hit.c) * std::log(p);
  }
  return sum;
}

// One reciprocal. Four hits are dotted before any division so those
// divisions can stay in flight together.
inline void add_newton_hit(
  double p,
  double dp,
  double d2p,
  int count,
  double *g,
  double *hess
) {
  if (!(p > 1e-300)) {
    p = 1e-300;
  }
  const double inv = 1.0 / p;
  const double r = dp * inv;
  const double cd = static_cast<double>(count);
  *g += cd * r;
  *hess += cd * (d2p * inv - r * r);
}

void lane_score(
  const NewtonLane &lane,
  const SpectralRows &rows,
  const double *e,
  const double *de,
  const double *d2e,
  double *g_out,
  double *h_out
) {
  double g = 0.0;
  double hess = 0.0;
  int n = 0;
  const int n_hit = lane.n_hit;
  for (; n + 4 <= n_hit; n += 4) {
    double p[4];
    double dp[4];
    double d2p[4];
    for (int slot = 0; slot < 4; ++slot) {
      const NewtonHit &hit = lane.hits[n + slot];
      dot_eigen3(
        rows.v[hit.i],
        rows.u[hit.j],
        e,
        de,
        d2e,
        &p[slot],
        &dp[slot],
        &d2p[slot]
      );
    }
    for (int slot = 0; slot < 4; ++slot) {
      add_newton_hit(
        p[slot],
        dp[slot],
        d2p[slot],
        lane.hits[n + slot].c,
        &g,
        &hess
      );
    }
  }
  for (; n < n_hit; ++n) {
    const NewtonHit &hit = lane.hits[n];
    double p = 0.0;
    double dp = 0.0;
    double d2p = 0.0;
    dot_eigen3(rows.v[hit.i], rows.u[hit.j], e, de, d2e, &p, &dp, &d2p);
    add_newton_hit(p, dp, d2p, hit.c, &g, &hess);
  }
  *g_out = g;
  *h_out = hess;
}

void newton_batch(
  NewtonLane *lanes,
  const SpectralRows &rows,
  const EvoModel &model
) {
  double e[kNewtonBatch][EvoModel::kState];
  double de[kNewtonBatch][EvoModel::kState];
  double d2e[kNewtonBatch][EvoModel::kState];
  hoist_exponentials(lanes, model, e, de, d2e);
  for (int lane = 0; lane < kNewtonBatch; ++lane) {
    if (!lanes[lane].present) {
      continue;
    }
    lanes[lane].grid_ll = lane_loglik(lanes[lane], rows, e[lane]);
    lanes[lane].out = lanes[lane].grid_t;
  }
  for (int iter = 0; iter < 8; ++iter) {
    bool any = false;
    for (int lane = 0; lane < kNewtonBatch; ++lane) {
      any = any || lanes[lane].moving;
    }
    if (!any) {
      break;
    }
    hoist_exponentials(lanes, model, e, de, d2e);
    for (int lane = 0; lane < kNewtonBatch; ++lane) {
      if (!lanes[lane].moving) {
        continue;
      }
      double g = 0.0;
      double h = 0.0;
      lane_score(lanes[lane], rows, e[lane], de[lane], d2e[lane], &g, &h);
      if (!std::isfinite(g) || !std::isfinite(h) || !(h < -1e-8)) {
        lanes[lane].moving = false;
        continue;
      }
      double t2 = lanes[lane].t - g / h;
      if (t2 <= lanes[lane].t_lo) {
        t2 = 0.5 * (lanes[lane].t + lanes[lane].t_lo);
      }
      if (t2 >= lanes[lane].t_hi) {
        t2 = 0.5 * (lanes[lane].t + lanes[lane].t_hi);
      }
      if (!std::isfinite(t2)) {
        lanes[lane].moving = false;
        continue;
      }
      const double step = std::fabs(t2 - lanes[lane].t);
      lanes[lane].t = t2;
      if (step < 1e-8 * (1.0 + t2)) {
        lanes[lane].moving = false;
      }
    }
  }
  hoist_exponentials(lanes, model, e, de, d2e);
  for (int lane = 0; lane < kNewtonBatch; ++lane) {
    if (!lanes[lane].present) {
      continue;
    }
    const double ll = lane_loglik(lanes[lane], rows, e[lane]);
    if (ll + 1e-8 >= lanes[lane].grid_ll) {
      lanes[lane].out = lanes[lane].t;
    }
  }
}

void prepare_lane(
  NewtonLane &lane,
  const EvoModel &model,
  int best
) {
  lane.present = true;
  lane.moving = true;
  lane.grid_t = model.t_grid[static_cast<std::size_t>(best)];
  lane.t = lane.grid_t;
  lane.t_lo = model.t_grid[static_cast<std::size_t>(best - 1)];
  lane.t_hi = model.t_grid[static_cast<std::size_t>(best + 1)];
  lane.out = lane.grid_t;
  collect_hits(lane);
}

double finish_distance(
  const PackedAaSequenceSet &pss,
  int i,
  int j,
  int n_ok,
  int best,
  int min_overlap,
  const EvoModel &model
) {
  const int n_grid = model.n_grid();
  if (n_ok < min_overlap || n_ok == 0 || best < 0 || best >= n_grid) {
    return model.t_max;
  }
  if (best == 0 || best == n_grid - 1) {
    return model.t_grid[static_cast<std::size_t>(best)];
  }
  NewtonLane lanes[kNewtonBatch];
  fill_counts(pss, i, j, model, lanes[0].counts);
  prepare_lane(lanes[0], model, best);
  const SpectralRows rows = spectral_rows(model);
  newton_batch(lanes, rows, model);
  return lanes[0].out;
}

void refine_evo_targets(
  const PackedAaSequenceSet &pss,
  int query_i,
  const int *target_idx,
  const int *n_ok,
  const int *best,
  int n_target,
  int min_overlap,
  const EvoModel &model,
  double *dist_out,
  int *n_ok_out
) {
  const int n_grid = model.n_grid();
  std::vector<int> interior;
  interior.reserve(static_cast<std::size_t>(n_target));
  for (int j = 0; j < n_target; ++j) {
    const int nok = n_ok[j];
    const int g = best[j];
    if (nok < min_overlap || nok == 0 || g < 0 || g >= n_grid) {
      n_ok_out[j] = 0;
      dist_out[j] = model.t_max;
      continue;
    }
    n_ok_out[j] = nok;
    if (g == 0 || g == n_grid - 1) {
      dist_out[j] = model.t_grid[static_cast<std::size_t>(g)];
      continue;
    }
    interior.push_back(j);
  }
  if (interior.empty()) {
    return;
  }
  const SpectralRows rows = spectral_rows(model);
  for (std::size_t base = 0; base < interior.size(); base += kNewtonBatch) {
    NewtonLane lanes[kNewtonBatch];
    const int n = std::min(
      kNewtonBatch,
      static_cast<int>(interior.size() - base)
    );
    for (int lane = 0; lane < n; ++lane) {
      const int j = interior[base + static_cast<std::size_t>(lane)];
      fill_counts(
        pss,
        query_i,
        target_idx[j],
        model,
        lanes[lane].counts
      );
      prepare_lane(lanes[lane], model, best[j]);
    }
    newton_batch(lanes, rows, model);
    for (int lane = 0; lane < n; ++lane) {
      const int j = interior[base + static_cast<std::size_t>(lane)];
      dist_out[j] = lanes[lane].out;
    }
  }
}

// Add one residue pair's contiguous grid row into a double accumulator.
inline void add_grid_row(const float *row, int n_grid, double *acc) {
  int g = 0;
#if OPTIMOTU_AA_KERNEL == 4
  for (; g + 8 <= n_grid; g += 8) {
    const __m256 v = _mm256_loadu_ps(row + g);
    _mm512_storeu_pd(
      acc + g,
      _mm512_add_pd(_mm512_loadu_pd(acc + g), _mm512_cvtps_pd(v))
    );
  }
#elif OPTIMOTU_AA_KERNEL == 3
  for (; g + 8 <= n_grid; g += 8) {
    const __m256 v = _mm256_loadu_ps(row + g);
    const __m128 lo = _mm256_castps256_ps128(v);
    const __m128 hi = _mm256_extractf128_ps(v, 1);
    _mm256_storeu_pd(
      acc + g,
      _mm256_add_pd(_mm256_loadu_pd(acc + g), _mm256_cvtps_pd(lo))
    );
    _mm256_storeu_pd(
      acc + g + 4,
      _mm256_add_pd(_mm256_loadu_pd(acc + g + 4), _mm256_cvtps_pd(hi))
    );
  }
#elif OPTIMOTU_AA_KERNEL == 2 && defined(__ARM_NEON)
  for (; g + 2 <= n_grid; g += 2) {
    const float64x2_t d = vcvt_f64_f32(vld1_f32(row + g));
    vst1q_f64(acc + g, vaddq_f64(vld1q_f64(acc + g), d));
  }
#elif OPTIMOTU_AA_KERNEL == 2
  for (; g + 2 <= n_grid; g += 2) {
    __m128i bits = _mm_setzero_si128();
    std::memcpy(&bits, row + g, sizeof(float) * 2);
    const __m128d d = _mm_cvtps_pd(_mm_castsi128_ps(bits));
    _mm_storeu_pd(acc + g, _mm_add_pd(_mm_loadu_pd(acc + g), d));
  }
#endif
  for (; g < n_grid; ++g) {
    acc[g] += static_cast<double>(row[g]);
  }
}

struct GridScan {
  int n_ok = 0;
  int best = 0;
  double best_ll = -std::numeric_limits<double>::infinity();
};

GridScan scan_pair_scalar(
  const PackedAaSequenceSet &pss,
  int i,
  int j,
  const EvoModel &model
) {
  GridScan out;
  const uint8_t *a = pss.row(i);
  const uint8_t *b = pss.row(j);
  const int lo = overlap_lo(pss, i, j);
  const int hi = overlap_hi(pss, i, j);
  const int n_grid = model.n_grid();
  std::vector<double> sum(static_cast<std::size_t>(n_grid), 0.0);
  for (int k = lo; k < hi; ++k) {
    const uint8_t ca = a[k];
    const uint8_t cb = b[k];
    if (!evo_site_ok(model, ca, cb)) {
      continue;
    }
    ++out.n_ok;
    add_grid_row(
      model.logp_row(ca, cb),
      n_grid,
      sum.data()
    );
  }
  for (int g = 0; g < n_grid; ++g) {
    const double ll = sum[static_cast<std::size_t>(g)];
    if (ll > out.best_ll) {
      out.best_ll = ll;
      out.best = g;
    }
  }
  return out;
}

}  // namespace

std::tuple<bool, int, double> PackedAaSequenceSet::score_evo_pair(
  int i,
  int j,
  int min_overlap,
  bool ignore_gaps,
  const EvoModel &model
) const {
  // Gaps and nonstandard residues have no state, so they are missing
  // data whether or not ignore_gaps is set.
  (void)ignore_gaps;
  const GridScan scan = scan_pair_scalar(*this, i, j, model);
  const double d = finish_distance(
    *this, i, j, scan.n_ok, scan.best, min_overlap, model
  );
  if (scan.n_ok < min_overlap || scan.n_ok == 0) {
    return {false, 0, model.t_max};
  }
  return {true, scan.n_ok, d};
}

double PackedAaSequenceSet::dist_evo(
  int i,
  int j,
  int min_overlap,
  bool ignore_gaps,
  const EvoModel &model
) const {
  return std::get<2>(score_evo_pair(i, j, min_overlap, ignore_gaps, model));
}

void PackedAaSequenceSet::score_query_evo(
  int query_i,
  const uint8_t *pos_major,
  const int *target_idx,
  int n_target,
  int n_cols,
  int min_overlap,
  bool ignore_gaps,
  const EvoModel &model,
  double *dist_out,
  int *n_ok_out
) const {
  (void)ignore_gaps;
  if (n_target <= 0) {
    return;
  }
  const uint8_t *q = row(query_i);
  const int q0 = start[static_cast<std::size_t>(query_i)];
  const int q1 = end[static_cast<std::size_t>(query_i)];
  bool covered = true;
  for (int j = 0; j < n_target; ++j) {
    const int tj = target_idx[j];
    if (start[static_cast<std::size_t>(tj)] > q0 ||
        end[static_cast<std::size_t>(tj)] < q1) {
      covered = false;
      break;
    }
  }

  const int n_grid = model.n_grid();
  std::vector<int> n_ok(static_cast<std::size_t>(n_target), 0);
  std::vector<int> best(static_cast<std::size_t>(n_target), 0);
  if (covered) {
    std::vector<double> acc(
      static_cast<std::size_t>(n_target) * static_cast<std::size_t>(n_grid),
      0.0
    );
    for (int k = q0; k < q1; ++k) {
      const uint8_t a = q[k];
      if (a >= EvoModel::kCode || model.is_state[a] == 0) {
        continue;
      }
      const float *page = model.logp_row(a, 0);
      const uint8_t *col = pos_major +
        static_cast<std::size_t>(k) * static_cast<std::size_t>(n_cols);
      for (int j = 0; j < n_target; ++j) {
        const uint8_t b = col[j];
        if (b >= EvoModel::kCode || model.is_state[b] == 0) {
          continue;
        }
        ++n_ok[static_cast<std::size_t>(j)];
        add_grid_row(
          page + static_cast<std::size_t>(b) * static_cast<std::size_t>(n_grid),
          n_grid,
          acc.data() +
            static_cast<std::size_t>(j) * static_cast<std::size_t>(n_grid)
        );
      }
    }
    for (int j = 0; j < n_target; ++j) {
      const double *row = acc.data() +
        static_cast<std::size_t>(j) * static_cast<std::size_t>(n_grid);
      double best_ll = -std::numeric_limits<double>::infinity();
      for (int g = 0; g < n_grid; ++g) {
        if (row[g] > best_ll) {
          best_ll = row[g];
          best[static_cast<std::size_t>(j)] = g;
        }
      }
    }
  } else {
    for (int j = 0; j < n_target; ++j) {
      const GridScan scan = scan_pair_scalar(
        *this,
        query_i,
        target_idx[j],
        model
      );
      n_ok[static_cast<std::size_t>(j)] = scan.n_ok;
      best[static_cast<std::size_t>(j)] = scan.best;
    }
  }
  refine_evo_targets(
    *this,
    query_i,
    target_idx,
    n_ok.data(),
    best.data(),
    n_target,
    min_overlap,
    model,
    dist_out,
    n_ok_out
  );
}
