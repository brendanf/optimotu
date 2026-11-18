// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_HAMMING_KERNEL_H
#define OPTIMOTU_HAMMING_KERNEL_H

// Compile-time Hamming stripe kernel. No target attributes.
// 1 = scalar, 2 = AVX2, 3 = AVX-512 VPOPCNTDQ, 4 = NEON.
#ifndef OPTIMOTU_HAMMING_KERNEL
#if defined(__AVX512F__) && defined(__AVX512DQ__) && defined(__AVX512VL__) && \
    defined(__AVX512VPOPCNTDQ__)
#define OPTIMOTU_HAMMING_KERNEL 3
#elif defined(__AVX2__)
#define OPTIMOTU_HAMMING_KERNEL 2
#elif defined(__ARM_NEON)
#define OPTIMOTU_HAMMING_KERNEL 4
#else
#define OPTIMOTU_HAMMING_KERNEL 1
#endif
#endif

#if OPTIMOTU_HAMMING_KERNEL == 3
#if !defined(__AVX512F__) || !defined(__AVX512VPOPCNTDQ__)
#error "OPTIMOTU_HAMMING_KERNEL=3 requires AVX-512F and VPOPCNTDQ"
#endif
#elif OPTIMOTU_HAMMING_KERNEL == 2
#if !defined(__AVX2__)
#error "OPTIMOTU_HAMMING_KERNEL=2 requires AVX2"
#endif
#elif OPTIMOTU_HAMMING_KERNEL == 4
#if !defined(__ARM_NEON)
#error "OPTIMOTU_HAMMING_KERNEL=4 requires ARM NEON"
#endif
#elif OPTIMOTU_HAMMING_KERNEL != 1
#error "OPTIMOTU_HAMMING_KERNEL must be 1, 2, 3, or 4"
#endif

// pdistB2 one-versus-many stripe: AVX2 / AVX-512 only. Scalar and NEON
// lose to pairwise pdistB2 on the CPUs we have timed.
#if OPTIMOTU_HAMMING_KERNEL == 2 || OPTIMOTU_HAMMING_KERNEL == 3
#define OPTIMOTU_HAMMING_PDISTB2_STRIPE 1
#else
#define OPTIMOTU_HAMMING_PDISTB2_STRIPE 0
#endif

// Word-major stripe for this ignore_gaps flag. Scalar (kernel 1, the
// Windows -msse2 build) stays on pairwise pdist / pdistB2. That is the
// path which matched hclust; the scalar stripe did not. AVX2, AVX-512,
// and NEON keep the stripe.
inline bool optimotu_hamming_use_stripe(bool ignore_gaps) {
#if OPTIMOTU_HAMMING_KERNEL == 1
  (void)ignore_gaps;
  return false;
#else
  return ignore_gaps || OPTIMOTU_HAMMING_PDISTB2_STRIPE;
#endif
}

#endif
