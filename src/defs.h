// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-FileContributor: Panu Somervuo
// SPDX-License-Identifier: MIT

#ifndef _DEFS_
#define _DEFS_

#include <stdint.h>

#define NUCLEOTIDES_IN_WORD 16

int nucleotide2binary(const char *s, const int n, uint64_t *b, uint64_t *m, int *start, int *end);

// Original from ProtaxA
// *a: first 4-bit encoded sequence
// *ma: first 1-bit encoded mask; 1 = base is valid
// *b: second 4-bit encoded sequence
// *mb: second 1-bit encoded mask; 1 = base is valid
// *start: start position to calculate the distance
// *end: end position to calculate the distance
// *min_len: minimum length of overlap between the two sequences
// *num_ok: number of bases that are valid in both sequences
// *num_matches: number of matches between the two sequences
// *dist: distance between the two sequences
// invalid characters (i.e., Ns, gaps, etc.) in either sequence are not
// counted as mismatches or as part of the alignment length
int pdistB(const uint64_t *a, const uint64_t *ma,
           const uint64_t *b, const uint64_t *mb,
           const int start, const int end, const int min_len,
           int *num_ok, int *num_matches, double *dist);

// As pdistB, but _internal_ invalid characters are treates as gaps;
// i.e. if a base is invalid in one of the sequences, it is counted as a
// mismatch, and also as part of the alignment length.
// However, if _both_ bases are invalid, the position is not counted as
// a mismatch or as part of the alignment length.
// Assuming that start and end give the positions of the first and last
// positions where both sequences are valid, then only "internal" gaps
// are counted.
int pdistB2(const uint64_t *a, const uint64_t *ma,
            const uint64_t *b, const uint64_t *mb,
            const int start, const int end, const int min_len,
            int *num_ok, int *num_matches, double *dist);

// As pdistB, but additionally calculates information about the internal gaps:
// *num_ins: number of insertions in the alignment (valid b, invalid a)
// *num_del: number of deletions in the alignment (valid a, invalid b)
// *max_ins: maximum length of an insertion
// *max_del: maximum length of a deletion
// Note that in pdistB, internal gaps do _not_ count in the calculation of
// the distance.
int pdistB_gap(const uint64_t *a, const uint64_t *ma,
               const uint64_t *b, const uint64_t *mb,
               const int start, const int end, const int min_len,
               int *num_ok, int *num_matches, double *dist,
               int *num_ins, int *num_del,
               int *max_ins, int *max_del);

// As pdistB2, but additionally calculates information about the internal gaps:
// *num_ins: number of insertions in the alignment (valid b, invalid a)
// *num_del: number of deletions in the alignment (valid a, invalid b)
// *max_ins: maximum length of an insertion
// *max_del: maximum length of a deletion
int pdistB2_gap(const uint64_t *a, const uint64_t *ma,
                const uint64_t *b, const uint64_t *mb,
                const int start, const int end, const int min_len,
                int *num_ok, int *num_matches, double *dist,
                int *num_ins, int *num_del,
                int *max_ins, int *max_del);

// As pdistB, but also calculates the CIGAR string of the alignment.
// *cigar: 0-terminated CIGAR string of the alignment. Should be preallocated
// to at least the alignment width
// Like pdistB_gap, this function is strange because it returns information
// about the gaps, but the gaps are ignored in calculating the distance.
// int pdistB_cigar(const uint64_t *a, const uint64_t *ma,
//                  const uint64_t *b, const uint64_t *mb,
//                  const int start, const int end, const int min_len,
//                  int *num_ok, int *num_matches, double *dist,
//                  char *cigar);

// As pdistB2, but also calculates the CIGAR string of the alignment.
// int pdistB2_cigar(const uint64_t *a, const uint64_t *ma,
//                   const uint64_t *b, const uint64_t *mb,
//                   const int start, const int end, const int min_len,
//                   int *num_ok, int *num_matches, double *dist,
//                   char *cigar);
#endif
