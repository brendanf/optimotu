#include "defs.h"

int nucleotide2binary(const char *s, const int n, long unsigned int *b, long unsigned int *m, int *start, int *end) {
  long unsigned int a, am;
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
    m[j] = am;
  }

  return (0);
}
#pragma GCC target ("sse4.2")
double pdistB(const long unsigned int *a, const long unsigned int *ma,
              const long unsigned int *b, const long unsigned int *mb,
              const int start, const int end, const int min_len)
{
  int i, num_ok, num_matches, nstart, nend;
  long unsigned int f;

  num_ok=0;
  num_matches=0;
  nstart=start/(NUCLEOTIDES_IN_WORD*4);
  nend=end/(NUCLEOTIDES_IN_WORD*4);
  if (nend * NUCLEOTIDES_IN_WORD*4 < end) nend++;

  for (i=nstart; i<nend; i++) {
    num_ok += __builtin_popcountl(ma[i] & mb[i]);
  }

  nstart=start/NUCLEOTIDES_IN_WORD;
  nend=end/NUCLEOTIDES_IN_WORD;
  if (nend * NUCLEOTIDES_IN_WORD < end) nend++;

  for (i=nstart; i<nend; i++) {
    num_matches += __builtin_popcountl(a[i] & b[i]);
  }

  if (num_ok >= min_len)
    return (1.0 - (double) num_matches / num_ok);
  else
    return (1.0);
}

#pragma GCC target ("sse4.2")
double pdistB2(const long unsigned int *a, const long unsigned int *ma,
              const long unsigned int *b, const long unsigned int *mb,
              const int start, const int end, const int min_len)
{
  int i, num_ok, num_matches, nstart, nend;
  long unsigned int mask;
  nstart=start/(NUCLEOTIDES_IN_WORD*4);
  nend=end/(NUCLEOTIDES_IN_WORD*4);
  if (nend * NUCLEOTIDES_IN_WORD*4 < end) nend++;

  if (nstart + 1 == nend) {
    mask = 0xffffffffffffffff << (start % (NUCLEOTIDES_IN_WORD*4));
    mask &= ~(0xffffffffffffffff << (end % (NUCLEOTIDES_IN_WORD*4)));
    num_ok = __builtin_popcountl((ma[nstart] | mb[nstart]) & mask);
  } else {
    mask = 0xffffffffffffffff << (start % (NUCLEOTIDES_IN_WORD*4));
    num_ok = __builtin_popcountl((ma[nstart] | mb[nstart]) & mask);

    for (i=nstart + 1; i<nend - 1; i++) {
      num_ok += __builtin_popcountl(ma[i] | mb[i]);
    }

    mask = ~(0xffffffffffffffff << (end % (NUCLEOTIDES_IN_WORD*4)));
    num_ok += __builtin_popcount((ma[nend-1] | mb[nend-1]) & mask);
  }

  if (num_ok == 0) return 1.0;

  nstart=start/NUCLEOTIDES_IN_WORD;
  nend=end/NUCLEOTIDES_IN_WORD;
  if (nend * NUCLEOTIDES_IN_WORD < end) nend++;

  num_matches=0;
  for (i=nstart; i<nend; i++) {
    num_matches += __builtin_popcountl(a[i] & b[i]);
  }

  if (num_ok >= min_len)
    return (1.0 - (double) num_matches / num_ok);
  else
    return (1.0);
}
