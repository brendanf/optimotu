#ifndef _DEFS_
#define _DEFS_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define NUCLEOTIDES_IN_WORD 16
#define MAXLINE 1024

typedef struct {
  int num_seqs, alen, ulen, mulen;
  char **id;
  long unsigned int **b, **m;
} SequenceSetB;

/* routines_sequence.c */

int nucleotide2binary(const char *s, int n, long unsigned int *b, long unsigned int *m);
SequenceSetB *read_aligned_sequencesB(char *filename);
double pdistB(const long unsigned int *a, const long unsigned int *ma, const long unsigned int *b, const long unsigned int *mb, const int n, const int n2);
int compute_distances(SequenceSetB *a, long unsigned int *b, long unsigned int *m, double *pdistances);

#endif
