// Originally from ProtaxA https://github.com/psomervuo/protaxA/blob/main/c2.zip
// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-FileContributor: Panu Somervuo
// SPDX-License-Identifier: MIT

#include "defs.h"
#include <stdio.h>
#include <inttypes.h>

//#define DEBUG 1

#ifdef DEBUG
#define debugprintf(...) fprintf(stderr, __VA_ARGS__)
#else
#define debugprintf(...)
#endif

int nucleotide2binary(const char *s, const int n, uint64_t *b, uint64_t *m, int *start, int *end) {
  uint64_t a, am;
  int i, // index in input sequence
      j, // index of output word
      k, // index within output word
      l, // index in output sequence
      n2, // number of words
      n_remaining, // number of bases in last incomplete word
      found_start = 0;

  /* sequence content, 4 bits for one character */

  n2 = n / NUCLEOTIDES_IN_WORD;
  i=0;
  for (j=0; j<n2; j++) {
    a = 0;
    for (k=0; k<NUCLEOTIDES_IN_WORD; k++) {
      a <<= 4;
      while ((s[i] >= 'a') && (s[i] <= 'z')) i++;
      if (s[i] == 'A') {a += 1;}
      else if (s[i] == 'C') {a += 2;}
      else if (s[i] == 'G') {a += 4;}
      else if (s[i] == 'T') {a += 8;}
      i++;
    }
    b[j] = a;
  }

  n_remaining = n - n2*NUCLEOTIDES_IN_WORD;
  if (n_remaining) {
    a = 0;
    for (k=0; k<n_remaining; k++) {
      a <<= 4;
      while ((s[i] >= 'a') && (s[i] <= 'z')) i++;
      if (s[i] == 'A') {a += 1;}
      else if (s[i] == 'C') {a += 2;}
      else if (s[i] == 'G') {a += 4;}
      else if (s[i] == 'T') {a += 8;}
      i++;
    }
    a <<= (4 * (NUCLEOTIDES_IN_WORD - n_remaining));
    b[j] = a;
  }

  /* mask, 1 bit for character: 1 ok, 0 not */

  n2 = n / (NUCLEOTIDES_IN_WORD*4);
  i=0;
  l=0;
  for (j=0; j<n2; j++) {
    am = 0;
    for (k=0; k<NUCLEOTIDES_IN_WORD*4; k++) {
      am <<= 1;
      while ((s[i] >= 'a') && (s[i] <= 'z')) {
        i++;
      }
      if ((s[i] == 'A') || (s[i] == 'C') || (s[i] == 'G') || (s[i] == 'T')) {
        if (!found_start) {
          found_start = 1;
          *start=l;
        }
        *end=l+1;
        am += 1;
      }
      i++;
      l++;
    }
    m[j] = am;
  }

  n_remaining = n - n2*NUCLEOTIDES_IN_WORD*4;
  if (n_remaining) {
    am = 0;
    for (k=0; k<n_remaining; k++) {
      am <<= 1;
      while ((s[i] >= 'a') && (s[i] <= 'z')) i++;
      if ((s[i] == 'A') || (s[i] == 'C') || (s[i] == 'G') || (s[i] == 'T')) {
        if (!found_start) {
          found_start = 1;
          *start=l;
        }
        *end=l+1;
        am += 1;
      }
      i++;
      l++;
    }
    am <<= (NUCLEOTIDES_IN_WORD*4 - n_remaining);
    m[j] = am;
  }

  return (0);
}

int pdistB(const uint64_t *a, const uint64_t *ma,
              const uint64_t *b, const uint64_t *mb,
              const int start, const int end, const int min_len,
              int *num_ok, int *num_matches, double *dist)
{
  int i, nstart, nend;

  debugprintf("entering pdistB:\n");
  debugprintf("a: ");
  for (i=0; i * NUCLEOTIDES_IN_WORD < end; i++) {
    debugprintf(" %016" PRIx64, a[i]);
  }
  debugprintf("\nma:");
  for (i=0; i * NUCLEOTIDES_IN_WORD*4 < end; i++) {
    debugprintf(" %016" PRIx64, ma[i]);
  }
  debugprintf("\nb: ");
  for (i=0; i * NUCLEOTIDES_IN_WORD < end; i++) {
    debugprintf(" %016" PRIx64, b[i]);
  }
  debugprintf("\nmb:");
  for (i=0; i * NUCLEOTIDES_IN_WORD*4 < end; i++) {
    debugprintf(" %016" PRIx64, mb[i]);
  }
  debugprintf("\n");

  *num_ok=0;
  *num_matches=0;
  nstart=start/(NUCLEOTIDES_IN_WORD*4);
  nend=end/(NUCLEOTIDES_IN_WORD*4);
  if (nend * NUCLEOTIDES_IN_WORD*4 < end) nend++;

  for (i=nstart; i<nend; i++) {
    *num_ok += __builtin_popcountll(ma[i] & mb[i]);
    debugprintf("i=%d, num_ok=%d\n", i, *num_ok);
  }

  if (*num_ok == 0 || *num_ok < min_len) {
    *dist = 1.0;
    return 1;
  }

  nstart=start/NUCLEOTIDES_IN_WORD;
  nend=end/NUCLEOTIDES_IN_WORD;
  if (nend * NUCLEOTIDES_IN_WORD < end) nend++;

  for (i=nstart; i<nend; i++) {
    *num_matches += __builtin_popcountll(a[i] & b[i]);
    debugprintf("i=%d, num_matches=%d\n", i, *num_matches);
  }

  *dist = 1.0 - (double) *num_matches / *num_ok;
  return 0;
}

int pdistB_gap(const uint64_t *a, const uint64_t *ma,
               const uint64_t *b, const uint64_t *mb,
               const int start, const int end, const int min_len,
               int *num_ok, int *num_matches, double *dist,
               int *num_ins, int *num_del,
               int *max_ins, int *max_del)
{
  int i, nstart, nend;

  debugprintf("entering pdistB_gap:\n");
  debugprintf("a: ");
  for (i=0; i * NUCLEOTIDES_IN_WORD < end; i++) {
    debugprintf(" %016" PRIx64, a[i]);
  }
  debugprintf("\nma:");
  for (i=0; i * NUCLEOTIDES_IN_WORD*4 < end; i++) {
    debugprintf(" %016" PRIx64, ma[i]);
  }
  debugprintf("\nb: ");
  for (i=0; i * NUCLEOTIDES_IN_WORD < end; i++) {
    debugprintf(" %016" PRIx64, b[i]);
  }
  debugprintf("\nmb:");
  for (i=0; i * NUCLEOTIDES_IN_WORD*4 < end; i++) {
    debugprintf(" %016" PRIx64, mb[i]);
  }
  debugprintf("\n");

  *num_ok=0;
  *num_matches=0;
  *num_ins=0;
  *num_del=0;
  *max_ins=0;
  *max_del=0;
  nstart=start/(NUCLEOTIDES_IN_WORD*4);
  nend=end/(NUCLEOTIDES_IN_WORD*4);
  if (nend * NUCLEOTIDES_IN_WORD*4 < end) nend++;

  int trail_ins = 0;
  int trail_del = 0;
  for (i=nstart; i<nend; i++) {
    unsigned long long end_mask = -1;
    if (i == nstart) {
      end_mask >>= (start % (NUCLEOTIDES_IN_WORD*4));
    }
    if (i == nend - 1) {
      end_mask &= ~((uint64_t)(-1) >> (end % (NUCLEOTIDES_IN_WORD*4)));
    }
    *num_ok += __builtin_popcountll((ma[i] & mb[i]) & end_mask);
    unsigned long long ins = mb[i] & ~ma[i] & end_mask;
    if (ins) {
      *num_ins += __builtin_popcountll(ins);
      int my_max_ins = 1;
      unsigned long long ins_mask = ins;
      while (ins_mask &= (ins_mask >> 1)) my_max_ins++;
      if (trail_ins && (ins_mask & 0x8000000000000000)) {
        trail_ins += __builtin_clrsbll(ins_mask) + 1;
      }
      my_max_ins = my_max_ins > trail_ins ? my_max_ins : trail_ins;
      *max_ins = *max_ins > my_max_ins ? *max_ins : my_max_ins;
      if (ins_mask != (unsigned long long)(-1)) {
        trail_ins = __builtin_ctzll(~ins_mask);
      }
    }
    unsigned long long del = ma[i] & ~mb[i] & end_mask;
    if (del) {
      *num_del += __builtin_popcountll(del);
      int my_max_del = 1;
      unsigned long long del_mask = del;
      while (del_mask &= (del_mask >> 1)) my_max_del++;
      if (trail_del && (del_mask & 0x8000000000000000)) {
        trail_del += __builtin_clrsbll(del_mask) + 1;
      }
      my_max_del = my_max_del > trail_del ? my_max_del : trail_del;
      *max_del = *max_del > my_max_del ? *max_del : my_max_del;
      if (del_mask != (unsigned long long)(-1)) {
        trail_del = __builtin_ctzll(~del_mask);
      }
    }
    debugprintf("i=%d, num_ok=%d, num_ins=%d, num_del=%d, max_ins=%d, max_del=%d, trail_ins=%d, trail_del=%d\n", i, *num_ok, *num_ins, *num_del, *max_ins, *max_del, trail_ins, trail_del);
  }

  if (*num_ok == 0 || *num_ok < min_len) {
    *dist = 1.0;
    return 1;
  }

  nstart=start/NUCLEOTIDES_IN_WORD;
  nend=end/NUCLEOTIDES_IN_WORD;
  if (nend * NUCLEOTIDES_IN_WORD < end) nend++;

  for (i=nstart; i<nend; i++) {
    *num_matches += __builtin_popcountll(a[i] & b[i]);
    debugprintf("i=%d, num_matches=%d\n", i, *num_matches);
  }

  *dist = 1.0 - (double) *num_matches / *num_ok;
  return 0;
}

int pdistB2(const uint64_t *a, const uint64_t *ma,
            const uint64_t *b, const uint64_t *mb,
            const int start, const int end, const int min_len,
            int *num_ok, int *num_matches, double *dist)
{
  int i, nstart, nend;
  uint64_t mask;

  debugprintf("entering pdistB2:\n");
  debugprintf("a: ");
  for (i=0; i * NUCLEOTIDES_IN_WORD < end; i++) {
    debugprintf(" %016" PRIx64, a[i]);
  }
  debugprintf("\nma:");
    for (i=0; i * (NUCLEOTIDES_IN_WORD*4) < end; i++) {
    debugprintf(" %016" PRIx64, ma[i]);
  }
  debugprintf("\nb: ");
  for (i=0; i * NUCLEOTIDES_IN_WORD < end; i++) {
    debugprintf(" %016" PRIx64, b[i]);
  }
  debugprintf("\nmb:");
  for (i=0; i * (NUCLEOTIDES_IN_WORD*4) < end; i++) {
    debugprintf(" %016" PRIx64, mb[i]);
  }
  debugprintf("\n");

  *num_ok = 0;
  *num_matches=0;
  nstart=start/(NUCLEOTIDES_IN_WORD*4);
  nend=end/(NUCLEOTIDES_IN_WORD*4);
  if (nend * NUCLEOTIDES_IN_WORD*4 < end) nend++;
  if (nstart + 1 == nend) {
    debugprintf("single mask\n");
    mask = (uint64_t)(-1) >> (start % (NUCLEOTIDES_IN_WORD*4));
    debugprintf("mask=%016" PRIx64 "\n", mask);
    mask &= ~((uint64_t)(-1) >> (end % (NUCLEOTIDES_IN_WORD*4)));
    debugprintf("mask=%016" PRIx64 "\n", mask);
    debugprintf("ma | mb=%016" PRIx64 "\n", ma[nstart] | mb[nstart]);
    *num_ok = __builtin_popcountll((ma[nstart] | mb[nstart]) & mask);
    debugprintf("num_ok=%d\n", *num_ok);
  } else {
    debugprintf("multiple masks\n");
    mask = (uint64_t)(-1) >> (start % (NUCLEOTIDES_IN_WORD*4));
    debugprintf("mask=%016" PRIx64 "\n", mask);
    *num_ok = __builtin_popcountll((ma[nstart] | mb[nstart]) & mask);
    debugprintf("num_ok=%d\n", *num_ok);

    for (i=nstart + 1; i<nend - 1; i++) {
      debugprintf("i=%d\n", i);
      debugprintf("ma | mb=%016" PRIx64 "\n", ma[i] | mb[i]);
      *num_ok += __builtin_popcountll(ma[i] | mb[i]);
      debugprintf("num_ok=%d\n", *num_ok);
    }

    mask = ~((uint64_t)(-1) >> (end % (NUCLEOTIDES_IN_WORD*4)));
    debugprintf("mask=%016" PRIx64 "\n", mask);
    *num_ok += __builtin_popcountll((ma[nend-1] | mb[nend-1]) & mask);
    debugprintf("num_ok=%d\n", *num_ok);
  }

  if (*num_ok == 0 || *num_ok < min_len) {
    *dist = 1.0;
    return 1;
  }

  nstart=start/NUCLEOTIDES_IN_WORD;
  nend=end/NUCLEOTIDES_IN_WORD;
  if (nend * NUCLEOTIDES_IN_WORD < end) nend++;

  for (i=nstart; i<nend; i++) {
    debugprintf("i=%d, ", i);
    *num_matches += __builtin_popcountll(a[i] & b[i]);
    debugprintf("num_matches=%d\n", *num_matches);
  }

  *dist = 1.0 - (double) *num_matches / *num_ok;
  return 0;
}

int pdistB2_gap(const uint64_t *a, const uint64_t *ma,
  const uint64_t *b, const uint64_t *mb,
  const int start, const int end, const int min_len,
  int *num_ok, int *num_matches, double *dist,
  int *num_ins, int *num_del,
  int *max_ins, int *max_del)
{
int i, nstart, nend;
uint64_t mask;

debugprintf("entering pdistB2_gap:\n");
debugprintf("a: ");
for (i=0; i * NUCLEOTIDES_IN_WORD < end; i++) {
debugprintf(" %016" PRIx64, a[i]);
}
debugprintf("\nma:");
for (i=0; i * (NUCLEOTIDES_IN_WORD*4) < end; i++) {
debugprintf(" %016" PRIx64, ma[i]);
}
debugprintf("\nb: ");
for (i=0; i * NUCLEOTIDES_IN_WORD < end; i++) {
debugprintf(" %016" PRIx64, b[i]);
}
debugprintf("\nmb:");
for (i=0; i * (NUCLEOTIDES_IN_WORD*4) < end; i++) {
debugprintf(" %016" PRIx64, mb[i]);
}
debugprintf("\n");

*num_ok = 0;
*num_matches=0;
nstart=start/(NUCLEOTIDES_IN_WORD*4);
nend=end/(NUCLEOTIDES_IN_WORD*4);
if (nend * NUCLEOTIDES_IN_WORD*4 < end) nend++;
if (nstart + 1 == nend) {
debugprintf("single mask\n");
mask = (uint64_t)(-1) >> (start % (NUCLEOTIDES_IN_WORD*4));
debugprintf("mask=%016" PRIx64 "\n", mask);
mask &= ~((uint64_t)(-1) >> (end % (NUCLEOTIDES_IN_WORD*4)));
debugprintf("mask=%016" PRIx64 "\n", mask);
debugprintf("ma | mb=%016" PRIx64 "\n", ma[nstart] | mb[nstart]);
*num_ok = __builtin_popcountll((ma[nstart] | mb[nstart]) & mask);
debugprintf("num_ok=%d\n", *num_ok);
unsigned long long ins = mb[nstart] & ~ma[nstart] & mask;
debugprintf("ins=%016" PRIx64 "\n", ins);
if (ins) {
  *num_ins += __builtin_popcountll(ins);
  int my_max_ins = 1;
  unsigned long long ins_mask = ins;
  while (ins_mask &= (ins_mask >> 1)) my_max_ins++;
}
unsigned long long del = ma[nstart] & ~mb[nstart] & mask;
debugprintf("del=%016" PRIx64 "\n", del);
if (del) {
  *num_del += __builtin_popcountll(del);
  int my_max_del = 1;
  unsigned long long del_mask = del;
  while (del_mask &= (del_mask >> 1)) my_max_del++;
}
debugprintf("num_ins=%d, num_del=%d, max_ins=%d, max_del=%d\n", *num_ins, *num_del, *max_ins, *max_del);
} else {
debugprintf("multiple masks\n");
// First mask
mask = (uint64_t)(-1) >> (start % (NUCLEOTIDES_IN_WORD*4));
debugprintf("mask=%016" PRIx64 "\n", mask);
*num_ok = __builtin_popcountll((ma[nstart] | mb[nstart]) & mask);
debugprintf("num_ok=%d\n", *num_ok);

unsigned long long ins = mb[nstart] & ~ma[nstart] & mask;
debugprintf("ins=%016" PRIx64 "\n", ins);
int trail_ins = 0;
if (ins) {
  *num_ins += __builtin_popcountll(ins);
  *max_ins = 1;
  unsigned long long ins_mask = ins;
  while (ins_mask &= (ins_mask >> 1)) *max_ins++;
  trail_ins = __builtin_ctzll(~ins_mask);
}
unsigned long long del = ma[nstart] & ~mb[nstart] & mask;
debugprintf("del=%016" PRIx64 "\n", del);
int trail_del = 0;
if (del) {
  *num_del += __builtin_popcountll(del);
  *max_del = 1;
  unsigned long long del_mask = del;
  while (del_mask &= (del_mask >> 1)) *max_del++;
  trail_del = __builtin_ctzll(~del_mask);
}
debugprintf("num_ins=%d, num_del=%d, max_ins=%d, max_del=%d\n", *num_ins, *num_del, *max_ins, *max_del);

// Middle masks
for (i=nstart + 1; i<nend - 1; i++) {
debugprintf("i=%d\n", i);
debugprintf("ma | mb=%016" PRIx64 "\n", ma[i] | mb[i]);
*num_ok += __builtin_popcountll(ma[i] | mb[i]);
debugprintf("num_ok=%d\n", *num_ok);
ins = mb[i] & ~ma[i];
debugprintf("ins=%016" PRIx64 "\n", ins);
if (ins) {
  *num_ins += __builtin_popcountll(ins);
  int my_max_ins = 1;
  unsigned long long ins_mask = ins;
  while (ins_mask &= (ins_mask >> 1)) my_max_ins++;
  if (trail_ins && (ins_mask & 0x8000000000000000)) {
    trail_ins += __builtin_clrsbll(ins_mask) + 1;
    my_max_ins = my_max_ins > trail_ins ? my_max_ins : trail_ins;
  }
  *max_ins = *max_ins > my_max_ins ? *max_ins : my_max_ins;
  if (ins_mask != (unsigned long long)(-1)) {
    trail_ins = __builtin_ctzll(~ins_mask);
  }
}
del = ma[i] & ~mb[i];
debugprintf("del=%016" PRIx64 "\n", del);
if (del) {
  *num_del += __builtin_popcountll(del);
  int my_max_del = 1;
  unsigned long long del_mask = del;
  while (del_mask &= (del_mask >> 1)) my_max_del++;
  if (trail_del && (del_mask & 0x8000000000000000)) {
    trail_del += __builtin_clrsbll(del_mask) + 1;
    my_max_del = my_max_del > trail_del ? my_max_del : trail_del;
  }
  *max_del = *max_del > my_max_del ? *max_del : my_max_del;
  if (del_mask != (unsigned long long)(-1)) {
    trail_del = __builtin_ctzll(~del_mask);
  }
}
}

// Last mask
mask = ~((uint64_t)(-1) >> (end % (NUCLEOTIDES_IN_WORD*4)));
debugprintf("mask=%016" PRIx64 "\n", mask);
*num_ok += __builtin_popcountll((ma[nend-1] | mb[nend-1]) & mask);
debugprintf("num_ok=%d\n", *num_ok);
ins = mb[nend-1] & ~ma[nend-1] & mask;
debugprintf("ins=%016" PRIx64 "\n", ins);
if (ins) {
  *num_ins += __builtin_popcountll(ins);
  int my_max_ins = 1;
  unsigned long long ins_mask = ins;
  while (ins_mask &= (ins_mask >> 1)) my_max_ins++;
  if (trail_ins && (ins_mask & 0x8000000000000000)) {
    trail_ins += __builtin_clrsbll(ins_mask) + 1;
    my_max_ins = my_max_ins > trail_ins ? my_max_ins : trail_ins;
  }
  *max_ins = *max_ins > my_max_ins ? *max_ins : my_max_ins;
}
del = ma[nend-1] & ~mb[nend-1] & mask;
debugprintf("del=%016" PRIx64 "\n", del);
if (del) {
  *num_del += __builtin_popcountll(del);
  int my_max_del = 1;
  unsigned long long del_mask = del;
  while (del_mask &= (del_mask >> 1)) my_max_del++;
  if (trail_del && (del_mask & 0x8000000000000000)) {
    trail_del += __builtin_clrsbll(del_mask) + 1;
    my_max_del = my_max_del > trail_del ? my_max_del : trail_del;
  }
  *max_del = *max_del > my_max_del ? *max_del : my_max_del;
}
debugprintf("num_ins=%d, num_del=%d, max_ins=%d, max_del=%d\n", *num_ins, *num_del, *max_ins, *max_del);
}

if (*num_ok == 0 || *num_ok < min_len) {
*dist = 1.0;
return 1;
}

nstart=start/NUCLEOTIDES_IN_WORD;
nend=end/NUCLEOTIDES_IN_WORD;
if (nend * NUCLEOTIDES_IN_WORD < end) nend++;

for (i=nstart; i<nend; i++) {
debugprintf("i=%d, ", i);
*num_matches += __builtin_popcountll(a[i] & b[i]);
debugprintf("num_matches=%d\n", *num_matches);
}

*dist = 1.0 - (double) *num_matches / *num_ok;
return 0;
}
