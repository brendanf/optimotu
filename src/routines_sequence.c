#include <zlib.h>
#include <pthread.h>
#include "defs.h"

#define WORKER_IDLE 0
#define WORKER_RUN 1
#define WORKER_STOP -1

int insertion_free_length(const char * line) {
  int i = 0, l = 0;
  while (line[i] != '\0' && line[i] != '\n' && i < MAXLINE) {
    if ((line[i] < 'a') || (line[i] > 'z')) {
      l++;
    }
    i++;
  };
  return l;
}

int insertion_free_copy(const char * line, char * target, int *start, int *end, const int len) {
  int i = 0, j = 0, found_start = 0;
  *start = 0;
  *end = 0;
  while (line[i] != '\0' && line[i] != '\n' && j < len && i < MAXLINE) {
    if ((line[i] < 'a') || (line[i] > 'z')) {
      target[j] = line[i];
      j++;
      if ((line[i] == 'A') || (line[i] == 'C') || (line[i] == 'G') || (line[i] == 'T')) {
        if (!found_start) {
          *start = j;
          found_start = 1;
        }
        *end = j + 1;
      }
    }
    i++;
  }
  return j;
}

void scan_aligned_sequences(const char *filename, int *alen, int *num_seqs) {
  gzFile fp;
  char line[MAXLINE];
  int linecount, i;
  char *thisfunction = "scan_aligned_sequences";

  if ((fp = gzopen(filename,"r")) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot open '%s' for reading.\n",thisfunction,filename);
    perror(""); exit(-1);
  }
  for (i=0;i<MAXLINE;i++) line[i] = '\0';
  linecount=0;
  while (gzgets(fp, line, MAXLINE)) {
    linecount++;
    if (line[MAXLINE-2] != '\0') {
      fprintf(stderr,"ERROR (%s): line %d length in file '%s' exceeds MAXLINE %d.\n",thisfunction,linecount,filename,MAXLINE);
      exit(-1);
    }
    if (line[0] != '>') {
      fprintf(stderr,"ERROR (%s): line %d in file '%s' doesn't start with '>' but '%c'.\n",thisfunction,linecount,filename,line[0]);
      exit(-1);
    }
    linecount++;
    if (gzgets(fp, line, MAXLINE) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot read line %d from file '%s'.\n",thisfunction,linecount,filename);
      exit(-1);
    }
    if (line[MAXLINE-2] != '\0') {
      fprintf(stderr,"ERROR (%s): line %d length in file '%s' exceeds MAXLINE %d.\n",thisfunction,linecount,filename,MAXLINE);
      exit(-1);
    }

  }

  /* calculate sequence length from last sequence */
  *alen = insertion_free_length(line);
  *num_seqs = linecount/2;

  gzclose(fp);
}

SequenceSet *read_aligned_sequences(const char *filename, const int len, const int num_seqs) {
  gzFile fp;
  char line[MAXLINE], *token;
  int i, iflen;
  SequenceSet *s;
  char *thisfunction = "read_aligned_sequences";

  if ((s = (SequenceSet *) malloc (sizeof(SequenceSet))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc SequenceSet.\n",thisfunction);
    perror(""); exit(-1);
  }

  s->alen = len;
  s->num_seqs = num_seqs;

  if ((fp = gzopen(filename,"r")) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot open '%s' for reading.\n",thisfunction,filename);
    perror(""); exit(-1);
  }

  for (i=0;i<MAXLINE;i++) line[i] = '\0';

  if ((s->id = (char **) malloc(s->num_seqs * sizeof(char *))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d char ptr array.\n",thisfunction,s->num_seqs);
    perror("");exit(-1);
  }
  if ((s->seq = (char **) malloc(s->num_seqs * sizeof(char *))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d*%d char array.\n",thisfunction,s->num_seqs,len+1);
    perror("");exit(-1);
  }
  if ((s->start = (int *) malloc(s->num_seqs * sizeof(int))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d int array.\n",thisfunction,s->num_seqs,len+1);
    perror("");exit(-1);
  }
  if ((s->end = (int *) malloc(s->num_seqs * sizeof(int))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d int array.\n",thisfunction,s->num_seqs,len+1);
    perror("");exit(-1);
  }
  for (i=0; i<s->num_seqs; i++) {
    if ((s->seq[i] = (char *) malloc((len+1) * sizeof(char))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d*%d char array.\n",thisfunction,s->num_seqs,len+1);
      perror("");exit(-1);
    }
  }

  for (i=0;i<num_seqs; i++) {
    if (gzgets(fp, line, MAXLINE) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot read entry %d name (linecount %d) from file '%s'.\n",thisfunction,i,2*i+1,filename);
      perror("");exit(-1);
    }
    token = strtok(line," \t\n|");
    s->id[i] = strdup(token+1);

    if (gzgets(fp, line, MAXLINE) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot read entry %d sequence (linecount %d) from file '%s'.\n",thisfunction,i,2*i+2,filename);
      perror("");exit(-1);
    }

    iflen = insertion_free_copy(line, s->seq[i], s->start + i, s->end + i, len);
    if (iflen != s->alen) {
      fprintf(stderr,"ERROR (%s): sequence lengths differ from %d, line %d, file '%s'.\n",thisfunction,len,2*i+2,filename);
      perror("");exit(-1);
    }
    s->seq[i][len] = '\0';
  }

  gzclose(fp);

  return(s);
}

double pdist(const char *a, const char *b, const int start, const int end,
             const int min_len) {
  int i,mismatches=0,okpositions=0;

  for (i=start; i<end; i++) {
    /* if (((a[i] == 'A') || (a[i] == 'C') || (a[i] == 'G') || (a[i] == 'T')) || ((b[i] == 'A') || (b[i] == 'C') || (b[i] == 'G') || (b[i] == 'T'))) { */
    if (((a[i] == 'A') || (a[i] == 'C') || (a[i] == 'G') || (a[i] == 'T')) && ((b[i] == 'A') || (b[i] == 'C') || (b[i] == 'G') || (b[i] == 'T'))) {
      okpositions++;
      if (a[i] != b[i])
        mismatches++;
    }
  }
  if (okpositions >= min_len)
    return ((double) mismatches/okpositions);
  else
    return (1.0);
}


int compute_distances(const SequenceSet *a, const char *seq, const int start,
                      const int end, const int min_len, double *pdistances)
{
  int i, s, e;

  if (end - start < min_len) {
    for (i = 0; i < a->num_seqs; i++) pdistances[i] = 1.0;
    return 0;
  }

  for (i=0; i<a->num_seqs; i++) {
    s = (start > a->start[i]) ? start : a->start[i];
    e = (end < a->end[i]) ? end : a->end[i];
    if (e - s < min_len) {
      pdistances[i] = 1.0;
    } else {
      pdistances[i] = pdist(seq, a->seq[i], s, e, min_len);
    }
  }
  return (0);
}


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

SequenceSetB * allocate_SequenceSetB(const int len, const int num_seqs) {
  SequenceSetB * s;
  int i, ulen, mulen;
  const char* thisfunction = "allocate_SequenceSetB";

  if ((s = (SequenceSetB *) malloc (sizeof(SequenceSetB))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc SequenceSetB.\n",thisfunction);
    perror(""); exit(-1);
  }

  s->alen = len;
  ulen = len / NUCLEOTIDES_IN_WORD;
  if (len > ulen * NUCLEOTIDES_IN_WORD)
    ulen++;
  s->ulen = ulen;
  mulen = len / NUCLEOTIDES_IN_WORD / 4;
  if (len > mulen * NUCLEOTIDES_IN_WORD / 4)
    mulen++;
  s->mulen = mulen;
  s->num_seqs = num_seqs;

  if ((s->id = (char **) malloc(s->num_seqs * sizeof(char *))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d char ptr array.\n",thisfunction,s->num_seqs);
    perror("");exit(-1);
  }
  if ((s->b = (long unsigned int **) malloc(s->num_seqs * sizeof(long unsigned int *))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d*%d ul array.\n",thisfunction,s->num_seqs,ulen);
    perror("");exit(-1);
  }
  if ((s->m = (long unsigned int **) malloc(s->num_seqs * sizeof(long unsigned int *))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d*%d ul array.\n",thisfunction,s->num_seqs,ulen);
    perror("");exit(-1);
  }
  if ((s->start = (int *) malloc(s->num_seqs * sizeof(int))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d int array.\n",thisfunction,s->num_seqs);
    perror("");exit(-1);
  }
  if ((s->end = (int *) malloc(s->num_seqs * sizeof(int))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d int array.\n",thisfunction,s->num_seqs);
    perror("");exit(-1);
  }
  for (i=0; i<s->num_seqs; i++) {
    if ((s->b[i] = (long unsigned int *) malloc(ulen * sizeof(long unsigned int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d*%d ul array.\n",thisfunction,s->num_seqs,ulen);
      perror("");exit(-1);
    }
    if ((s->m[i] = (long unsigned int *) malloc(mulen * sizeof(long unsigned int))) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot malloc %d*%d ul array.\n",thisfunction,s->num_seqs,mulen);
      perror("");exit(-1);
    }
  }
  return s;
}

// allocates two SequenceSetB, such that their array members are contiguous
// arrays, meaning that values in the query sequence set can be reached by
// indexing past the "end" of the reference sequence set.
void allocate_SequenceSetB_pair(const int len, const int num_rseqs,
                                    const int num_qseqs,
                                    SequenceSetB ** r, SequenceSetB ** q) {
  int i, num_seqs, ulen, mulen;
  long unsigned int *b, *m;
  const char* thisfunction = "allocate_SequenceSetB";

  if ((*r = (SequenceSetB *) malloc (sizeof(SequenceSetB))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc SequenceSetB.\n",thisfunction);
    perror(""); exit(-1);
  }

  if ((*q = (SequenceSetB *) malloc (sizeof(SequenceSetB))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc SequenceSetB.\n",thisfunction);
    perror(""); exit(-1);
  }

  (*r)->alen = (*q)->alen = len;
  ulen = len / NUCLEOTIDES_IN_WORD;
  if (len > ulen * NUCLEOTIDES_IN_WORD)
    ulen++;
  (*r)->ulen = (*q)->ulen = ulen;
  mulen = len / NUCLEOTIDES_IN_WORD / 4;
  if (len > mulen * NUCLEOTIDES_IN_WORD / 4)
    mulen++;
  (*r)->mulen = (*q)->mulen = mulen;
  (*r)->num_seqs = num_rseqs;
  (*q)->num_seqs = num_qseqs;
  num_seqs = num_rseqs + num_qseqs;

  if (((*r)->id = (char **) malloc(num_seqs * sizeof(char *))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d char ptr array.\n",thisfunction,num_seqs);
    perror("");exit(-1);
  }
  (*q)->id = (*r)->id + num_rseqs;

  if (((*r)->b = (long unsigned int **) malloc(num_seqs * sizeof(long unsigned int *))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d*%d ul array.\n",thisfunction,num_seqs,ulen);
    perror("");exit(-1);
  }
  (*q)->b = (*r)->b + num_rseqs;

  if (((*r)->m = (long unsigned int **) malloc(num_seqs * sizeof(long unsigned int *))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d*%d ul array.\n",thisfunction,num_seqs,ulen);
    perror("");exit(-1);
  }
  (*q)->m = (*r)->m + num_rseqs;

  if (((*r)->start = (int *) malloc(num_seqs * sizeof(int))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d int array.\n",thisfunction,num_seqs);
    perror("");exit(-1);
  }
  (*q)->start = (*r)->start + num_rseqs;
  if (((*r)->end = (int *) malloc(num_seqs * sizeof(int))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d int array.\n",thisfunction,num_seqs);
    perror("");exit(-1);
  }
  (*q)->end = (*r)->end + num_rseqs;
  if (((*r)->b[0] = (long unsigned int *) malloc(num_seqs * ulen * sizeof(long unsigned int))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d*%d ul array.\n",thisfunction,num_seqs,ulen);
    perror("");exit(-1);
  }
  if (((*r)->m[0] = (long unsigned int *) malloc(num_seqs * mulen * sizeof(long unsigned int))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc %d*%d ul array.\n",thisfunction,num_seqs,mulen);
    perror("");exit(-1);
  }
  b = (*r)->b[0];
  m = (*r)->m[0];
  for (i=1; i<num_seqs; i++) {
    (*r)->b[i] = (b += ulen);
    (*r)->m[i] = (m += mulen);
  }
  return;

}

void read_allocated_SequenceSetB(SequenceSetB * s, const char *filename, const int len, const int num_seqs) {
  gzFile fp;
  char line[MAXLINE], *token;
  int i, iflen;
  char *thisfunction = "read_allocated_SequenceSetB";

  if ((fp = gzopen(filename,"r")) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot open '%s' for reading.\n",thisfunction,filename);
    perror(""); exit(-1);
  }

  for (i=0;i<MAXLINE;i++) line[i] = '\0';

  for (i=0;i<s->num_seqs; i++) {
    if (gzgets(fp, line, MAXLINE) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot read entry %d name (linecount %d) from file '%s'.\n",thisfunction,i,2*i+1,filename);
      perror("");exit(-1);
    }
    token = strtok(line," \t\n|");
    s->id[i] = strdup(token+1);

    if (gzgets(fp, line, MAXLINE) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot read entry %d sequence (linecount %d) from file '%s'.\n",thisfunction,i,2*i+2,filename);
      perror("");exit(-1);
    }

    iflen = insertion_free_length(line);
    if (iflen != len) {
      fprintf(stderr,"ERROR (%s): sequence length %d differs from %d, line %d, file '%s':\n%s\n",thisfunction,iflen,len,2*i+2,filename,line);
      exit(-1);
    }

    nucleotide2binary(line, len, s->b[i], s->m[i], s->start + i, s->end + i);
  }

  gzclose(fp);
}

SequenceSetB *read_aligned_sequencesB(const char *filename, const int len, const int num_seqs) {
  gzFile fp;
  char line[MAXLINE], *token;
  SequenceSetB *s;
  char *thisfunction = "read_aligned_sequencesB";

  s = allocate_SequenceSetB(len, num_seqs);
  read_allocated_SequenceSetB(s, filename, len, num_seqs);

  return(s);
}

void assert_length_match(const char *fname, const int len1, const int len2) {
  if (len1 != len2) {
    fprintf(
      stderr,
      "ERROR: sequence lengths in file %s do not match expected value (%d,%d).\n",
      fname,len1,len2
    );
    exit(1);
  }
}

void read_sequence_sets(InputOptions iopt, const char * rfile, const char * ifile,
                         SequenceSet **rseq, SequenceSet **iseq
) {
  int rlen, ilen;

  if (iopt.len > 0) {
    if (iopt.n_rseq > 0) {
      rlen = iopt.len;
    } else {
      scan_aligned_sequences(rfile, &rlen, &iopt.n_rseq);
      assert_length_match(rfile, iopt.len, rlen);
    }

    if (iopt.n_iseq > 0) {
      ilen = iopt.len;
    } else {
      scan_aligned_sequences(ifile, &ilen, &iopt.n_iseq);
      assert_length_match(ifile, iopt.len, ilen);
    }
  } else if (iopt.n_rseq > 0) {
    scan_aligned_sequences(ifile, &ilen, &iopt.n_iseq);
    rlen = ilen;
  } else {
    scan_aligned_sequences(rfile, &rlen, &iopt.n_rseq);
    if (iopt.n_iseq > 0) {
      ilen = rlen;
    } else {
      scan_aligned_sequences(ifile, &ilen, &iopt.n_iseq);
      assert_length_match(ifile, rlen, ilen);
    }
  }

  *rseq = read_aligned_sequences(rfile, rlen, iopt.n_rseq);
  *iseq = read_aligned_sequences(ifile, ilen, iopt.n_iseq);
}

void read_sequence_setsB(InputOptions iopt, const char * rfile, const char * ifile,
    SequenceSetB **rseq, SequenceSetB **iseq
) {
  int rlen, ilen;

  if (iopt.len > 0) {
    if (iopt.n_rseq > 0) {
      rlen = iopt.len;
    } else {
      scan_aligned_sequences(rfile, &rlen, &iopt.n_rseq);
      assert_length_match(rfile, iopt.len, rlen);
    }

    if (iopt.n_iseq > 0) {
      ilen = iopt.len;
    } else {
      scan_aligned_sequences(ifile, &ilen, &iopt.n_iseq);
      assert_length_match(ifile, iopt.len, ilen);
    }
  } else if (iopt.n_rseq > 0) {
    scan_aligned_sequences(ifile, &ilen, &iopt.n_iseq);
    rlen = ilen;
  } else {
    scan_aligned_sequences(rfile, &rlen, &iopt.n_rseq);
    if (iopt.n_iseq > 0) {
      ilen = rlen;
    } else {
      scan_aligned_sequences(ifile, &ilen, &iopt.n_iseq);
      assert_length_match(ifile, rlen, ilen);
    }
  }

  allocate_SequenceSetB_pair(rlen, iopt.n_rseq, iopt.n_iseq, rseq, iseq);
  read_allocated_SequenceSetB(*rseq, rfile, rlen, iopt.n_rseq);
  read_allocated_SequenceSetB(*iseq, ifile, ilen, iopt.n_iseq);
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
              const int n, const int n2, const int start, const int end)
{
  int i, num_ok, num_matches, nstart, nend;
  long unsigned int mask;
  fprintf(stderr, "start=%d end=%d ", start, end);
  nstart=start/(NUCLEOTIDES_IN_WORD*4);
  nend=end/(NUCLEOTIDES_IN_WORD*4);
  if (nend * NUCLEOTIDES_IN_WORD*4 < end) nend++;

  if (nstart + 1 == nend) {
    mask = 0xffffffffffffffff >> (start % (NUCLEOTIDES_IN_WORD*4));
    mask &= ~(0xffffffffffffffff >> (end % (NUCLEOTIDES_IN_WORD*4)));
    num_ok = __builtin_popcountl((ma[nstart] | mb[nstart]) & mask);
  } else {

    mask = 0xffffffffffffffff >> (start % (NUCLEOTIDES_IN_WORD*4));
    num_ok = __builtin_popcountl((ma[i] | mb[i]) & mask);

    for (i=nstart + 1; i<nend - 1; i++) {
      num_ok += __builtin_popcountl(ma[i] | mb[i]);
    }

    mask &= ~(0xffffffffffffffff >> (end % (NUCLEOTIDES_IN_WORD*4)));
    num_ok -= __builtin_popcount((ma[nend] | mb[nend]) & mask);
  }

  if (num_ok == 0) return 1.0;

  nstart=start/NUCLEOTIDES_IN_WORD;
  nend=NUCLEOTIDES_IN_WORD;
  if (nend * NUCLEOTIDES_IN_WORD < end) nend++;

  num_matches=0;
  for (i=nstart; i<nend; i++) {
    num_matches += __builtin_popcountl(a[i] & b[i]);
  }

  return (1.0 - (double) num_matches / num_ok);
}

int compute_distancesB_range(
    const SequenceSetB *a,
    const int a_from, const int a_to,
    const long unsigned int *b, const long unsigned int *m,
    const int start, const int end,
    const int min_len,
    double *pdistances
) {
  int i, s, e;
  if (end - start < min_len) {
    for (i = a_from; i < a_to; i++) pdistances[i] = 1.0;
    return 0;
  }

  for (i=a_from; i<a_to; i++) {
    s = (start > a->start[i]) ? start : a->start[i];
    e = (end < a->end[i]) ? end : a->end[i];
    if (e - s < min_len) {
      pdistances[i] = 1.0;
    } else {
      pdistances[i] = pdistB(b, m, a->b[i], a->m[i], s, e, min_len);
    }
  }
  return 0;
}

int compute_distancesB(
    const SequenceSetB *a,
    const long unsigned int *b, const long unsigned int *m,
    const int start, const int end,
    const int min_len,
    double *pdistances
) {
  compute_distancesB_range(a, 0, a->num_seqs, b, m, start, end, min_len, pdistances);
  return 0;
}



void * distancesB_worker_run(void * arg) {
  distancesB_worker *config = arg;

  while (1) {
    pthread_mutex_lock(&(config->mutex));

    while (config->worker_command == WORKER_IDLE) {
      pthread_cond_wait(&(config->condition), &(config->mutex));
    }

    pthread_mutex_unlock(&(config->mutex));

    if (config->worker_command == WORKER_STOP) break;

    compute_distancesB_range(
      config->a,
      config->a_from, config->a_to,
      config->b, config->m,
      config->start, config->end,
      config->min_len,
      config->pdistances
    );

    pthread_mutex_lock(&(config->mutex));
    if (config->worker_command == WORKER_STOP) {
      pthread_mutex_unlock(&(config->mutex));
      break;
    }
    config->worker_command = WORKER_IDLE;
    pthread_cond_signal(&(config->condition));
    pthread_mutex_unlock(&(config->mutex));
  }

  return NULL;
}

distancesB_worker * init_distancesB_worker_pool(
    int nthread, const SequenceSetB *a, const int min_len, double *pdistances
) {
  int i;
  pthread_t thread;
  const char* thisfunction = "init_distancesB_worker_pool";
  distancesB_worker * pool;
  if ((pool = malloc(nthread * sizeof(distancesB_worker))) == NULL) {
    fprintf(stderr, "ERROR: (%s) could not allocate pool of %d worker threads.\n",
            thisfunction, nthread);
    perror("");exit(-1);
  }
  for (i = 0; i < nthread; i++) {
    pool[i].a = a;
    pool[i].min_len = min_len;
    pool[i].pdistances = pdistances;
    pool[i].worker_command = WORKER_IDLE;
    pthread_mutex_init(&(pool[i].mutex), NULL);
    pthread_cond_init(&(pool[i].condition), NULL);
    pthread_create(&thread, NULL, distancesB_worker_run, &(pool[i]));
    pthread_detach(thread);
  }

  return pool;
}

int compute_distancesB_parallel(
    distancesB_worker * pool, int nthread,
    const long unsigned int *b, const long unsigned int *m,
    const int start, const int end
) {
  int i, a_to = 0;
  if (nthread < 1) nthread = 1;
  for (i = 0; i < nthread; i++) {
    pthread_mutex_lock(&(pool[i].mutex));
    pool[i].a_from = a_to;
    a_to = (pool[i].a->num_seqs * (i + 1)) / nthread;
    pool[i].a_to = a_to;
    // fprintf(stderr, "thread %d; %d to %d\n", i, pool[i].a_from, pool[i].a_to);
    pool[i].start = start;
    pool[i].end = end;
    pool[i].b = b;
    pool[i].m = m;
    pool[i].worker_command = WORKER_RUN;
    pthread_cond_signal(&(pool[i].condition));
    pthread_mutex_unlock(&(pool[i].mutex));
  }

  for (i = 0; i < nthread; i++) {
    pthread_mutex_lock(&(pool[i].mutex));
    while (pool[i].worker_command == WORKER_RUN) {
      pthread_cond_wait(&(pool[i].condition), &(pool[i].mutex));
    }
    pthread_mutex_unlock(&(pool[i].mutex));
  }
  return 0;
}

int stop_distancesB_worker_pool(distancesB_worker* pool, int nthread) {
  int i;

  for (int i = 0; i < nthread; i++) {
    pthread_mutex_lock(&(pool[i].mutex));
    pool[i].worker_command = WORKER_STOP;
    pthread_cond_signal(&(pool[i].condition));
    pthread_mutex_unlock(&(pool[i].mutex));
  }
  return 0;
}
