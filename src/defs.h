// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-FileContributor: Panu Somervuo
// SPDX-License-Identifier: MIT

#ifndef _DEFS_
#define _DEFS_

#define NUCLEOTIDES_IN_WORD 16

int nucleotide2binary(const char *s, const int n, long unsigned int *b, long unsigned int *m, int *start, int *end);
double pdistB(const long unsigned int *a, const long unsigned int *ma,
              const long unsigned int *b, const long unsigned int *mb,
              const int start, const int end, const int min_len);
double pdistB2(const long unsigned int *a, const long unsigned int *ma,
              const long unsigned int *b, const long unsigned int *mb,
              const int start, const int end, const int min_len);

#endif
