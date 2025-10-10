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
