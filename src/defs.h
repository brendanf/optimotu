// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-FileContributor: Panu Somervuo
// SPDX-License-Identifier: MIT

#ifndef _DEFS_
#define _DEFS_

#include <stdint.h>

#define NUCLEOTIDES_IN_WORD 16

int nucleotide2binary(const char *s, const int n, uint64_t *b, uint64_t *m, int *start, int *end);
double pdistB(const uint64_t *a, const uint64_t *ma,
              const uint64_t *b, const uint64_t *mb,
              const int start, const int end, const int min_len);
double pdistB2(const uint64_t *a, const uint64_t *ma,
              const uint64_t *b, const uint64_t *mb,
              const int start, const int end, const int min_len);

#endif
