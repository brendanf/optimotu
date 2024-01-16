#include "defs.h"

int nucleotide2binary(const char *s, int n, long unsigned int *b, long unsigned int *m) {
  long unsigned int a, am;
  int i,j,k, n2, n_remaining;

  /* sequence content, 4 bits for one character */

  n2 = n / NUCLEOTIDES_IN_WORD;
  i=0;
  for (j=0; j<n2; j++) {
    a = 0;
    for (k=0; k<NUCLEOTIDES_IN_WORD; k++) {
      a <<= 4;
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
  for (j=0; j<n2; j++) {
    am = 0;
    for (k=0; k<NUCLEOTIDES_IN_WORD*4; k++) {
      am <<= 1;
      if ((s[i] == 'A') || (s[i] == 'C') || (s[i] == 'G') || (s[i] == 'T'))
	am += 1;
      i++;
    }
    m[j] = am;
  }

  n_remaining = n - n2*NUCLEOTIDES_IN_WORD*4;
  if (n_remaining) {
    am = 0;
    for (k=0; k<n_remaining; k++) {
      am <<= 1;
      if ((s[i] == 'A') || (s[i] == 'C') || (s[i] == 'G') || (s[i] == 'T'))
	am += 1;
      i++;
    }
    m[j] = am;
  }

  return (0);
}


SequenceSetB *read_aligned_sequencesB(char *filename) {
  FILE *fp;
  char line[MAXLINE], *token;
  int len,linecount, i, ok, ulen, mulen;
  SequenceSetB *s;
  char *thisfunction = "read_aligned_sequences";

  if ((s = (SequenceSetB *) malloc (sizeof(SequenceSetB))) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot malloc SequenceSetB.\n",thisfunction);
    perror(""); exit(-1);
  }

  if ((fp = fopen(filename,"r")) == NULL) {
    fprintf(stderr,"ERROR (%s): cannot open '%s' for reading.\n",thisfunction,filename);
    perror(""); exit(-1);
  }

  for (i=0;i<MAXLINE;i++) line[i] = '\0';

  linecount=0;
  while (fgets(line, MAXLINE, fp)) {
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
    if (fgets(line, MAXLINE, fp) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot read line %d from file '%s'.\n",thisfunction,linecount,filename);
      exit(-1);
    }
    if (line[MAXLINE-2] != '\0') {
      fprintf(stderr,"ERROR (%s): line %d length in file '%s' exceeds MAXLINE %d.\n",thisfunction,linecount,filename,MAXLINE);
      exit(-1);
    }
  }

  /* calculate sequence length from last sequence */
  len=0;
  while (line[len] != '\n') len++;

  s->alen = len;
  ulen = len / NUCLEOTIDES_IN_WORD;
  if (len > ulen * NUCLEOTIDES_IN_WORD)
    ulen++;
  s->ulen = ulen;
  mulen = len / NUCLEOTIDES_IN_WORD / 4;
  if (len > mulen * NUCLEOTIDES_IN_WORD / 4)
    mulen++;
  s->mulen = mulen;
  s->num_seqs = linecount/2;

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

  rewind(fp);
  for (i=0;i<s->num_seqs; i++) {
    if (fgets(line, MAXLINE, fp) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot read entry %d name (linecount %d) from file '%s'.\n",thisfunction,i,2*i+1,filename);
      perror("");exit(-1);
    }
    token = strtok(line," \t\n");
    s->id[i] = strdup(token+1);

    if (fgets(line, MAXLINE, fp) == NULL) {
      fprintf(stderr,"ERROR (%s): cannot read entry %d sequence (linecount %d) from file '%s'.\n",thisfunction,i,2*i+2,filename);
      perror("");exit(-1);
    }
    if (strlen(line) != (len+1)) {
      fprintf(stderr,"ERROR (%s): sequence lenghts differ from %d, line %d, file '%s'.\n",thisfunction,len,2*i+2,filename);
      perror("");exit(-1);
    }

    nucleotide2binary(line, len, s->b[i], s->m[i]);
  }

  fclose(fp);

  return(s);
}

#pragma GCC target ("sse4.2")
double pdistB(
    const long unsigned int *a,
    const long unsigned int *ma,
    const long unsigned int *b,
    const long unsigned int *mb,
    const int n,
    const int n2
) {
  int i, num_ok, num_matches;
  long unsigned int f;

  num_ok=0;
  num_matches=0;

  for (i=0; i<n2; i++) {
    num_ok += __builtin_popcountl(ma[i] & mb[i]);
  }

  for (i=0; i<n; i++) {
    num_matches += __builtin_popcountl(a[i] & b[i]);
  }

  if (num_ok > 0)
    return (1.0 - (double) num_matches / num_ok);
  else
    return (1.0);
}

int compute_distances(SequenceSetB *a, long unsigned int *b, long unsigned int *m, double *pdistances)
{
  int i;

  for (i=0; i<a->num_seqs; i++) {
    pdistances[i] = pdistB(b, m, a->b[i], a->m[i], a->ulen, a->mulen);
  }
  return (0);
}
