// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_AA_KERNEL_H
#define OPTIMOTU_AA_KERNEL_H

// Compile-time amino-acid identity width. No target attributes.
// 1 = scalar, 2 = SSE2/NEON 16, 3 = AVX2 32, 4 = AVX-512BW 64.
#ifdef OPTIMOTU_AA_SCALAR_ONLY
#undef OPTIMOTU_AA_KERNEL
#define OPTIMOTU_AA_KERNEL 1
#endif

#ifndef OPTIMOTU_AA_KERNEL
#if defined(__AVX512BW__)
#define OPTIMOTU_AA_KERNEL 4
#elif defined(__AVX2__)
#define OPTIMOTU_AA_KERNEL 3
#elif defined(__SSE2__) || defined(__ARM_NEON)
#define OPTIMOTU_AA_KERNEL 2
#else
#define OPTIMOTU_AA_KERNEL 1
#endif
#endif

#if OPTIMOTU_AA_KERNEL == 4
#if !defined(__AVX512BW__)
#error "OPTIMOTU_AA_KERNEL=4 requires AVX-512BW"
#endif
#elif OPTIMOTU_AA_KERNEL == 3
#if !defined(__AVX2__)
#error "OPTIMOTU_AA_KERNEL=3 requires AVX2"
#endif
#elif OPTIMOTU_AA_KERNEL == 2
#if !defined(__SSE2__) && !defined(__ARM_NEON)
#error "OPTIMOTU_AA_KERNEL=2 requires SSE2 or ARM NEON"
#endif
#elif OPTIMOTU_AA_KERNEL != 1
#error "OPTIMOTU_AA_KERNEL must be 1, 2, 3, or 4"
#endif

#endif
