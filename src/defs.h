// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-FileContributor: Panu Somervuo
// SPDX-License-Identifier: MIT

#ifndef _DEFS_
#define _DEFS_

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#define NUCLEOTIDES_IN_WORD 16
#define MAXLINE 10240
#define UNKNAME "unk"

typedef struct {
  double rth;
  int len, n_rseq, n_iseq, min_len, nthread;
} InputOptions;

typedef struct {
  int nid; // node id; unique identifier of this node, should be ascending
  int pid; // parent id; unique identifier of this node's parent
  int level; // taxonomic rank, starts at 0 for root.
  char *name; // name of this taxon, not (?) including parents
  char *mxcode; // MX code; unique identifier in some database?
  int num_cnodes; //number of child nodes
  int size_cnodes; // size of the child node array (may include room to grow)
  int *cnode_index; // child node array
  int num_rseqs; // number of reference sequences
  int size_rseqs; // size of reference seq array (may include room to grow)
  int *rseq_index; // reference sequence array
  int isunk; // 1 if this node represents "unknown"
  int no_rseqs; // 1 if this node has no reference sequences
  double prior; // prior probability for this node
  double prob; //calculated probability for this node
  double sumcprob_no_rseqs;
  int ind1, ind2; // index of refseqs with best and second best match to query
  double dist1, dist2; // distances to first and second best match to query
} TaxonomyNode;

typedef struct {
  TaxonomyNode *node; // actual taxonomy nodes
  int *discovered; // number of newly discovered taxa at each rank
  char **rank_name; //names of ranks
  int max_level; // number of ranks in the taxonomy
  int num_nodes; // number of nodes actually used
  int size; // number of nodes allocated (may be room to grow)
} Taxonomy;

typedef struct {
  int num_seqs, alen;
  char **seq, **id;
  int *start, *end;
} SequenceSet;

typedef struct {
  int num_seqs, alen, ulen, mulen;
  char **id;
  long unsigned int **b, **m;
  int *start, *end;
} SequenceSetB;

typedef struct {
  int num_levels, dim;
  double **params;
} Model;

typedef struct{
  const SequenceSetB *a;
  double *pdistances;
  pthread_mutex_t mutex;
  pthread_cond_t condition;
  int worker_command;
  int min_len;
  int a_from, a_to, start, end;
  const long unsigned int *b, *m;
} distancesB_worker;

/* routines_taxonomy.c */

Taxonomy read_taxonomy(char *filename);
void grow_taxonomy(Taxonomy *t);
int add_rseq2taxonomy(char *filename, TaxonomyNode *node);
void add_rank_names(char* rankarg, Taxonomy *t);
int print_taxonomy(TaxonomyNode *node, int num_nodes);
int add_cnode(Taxonomy t, int nid, int cid);
void name_as_pseudo(Taxonomy *t, TaxonomyNode *n);
void name_as_unknown(Taxonomy *t, TaxonomyNode *n);
int add_taxon(Taxonomy *t, int pid, int isunk);
int add_rseq(Taxonomy t, int nid, int seqid);
// returns the taxon id of the "known" leaf node
int fill_pseudotree(Taxonomy *t, int nid, int qi);

/* routines_cli.c */

InputOptions get_input_options(const int argc, char * const * argv);
InputOptions get_input_options_custom(const int argc, char * const * argv, const char * options);

/* routines_sequence.c */

void scan_aligned_sequences(const char *filename, int *len, int *num_seqs);

SequenceSet *read_aligned_sequences(const char *filename, const int len, const int num_seqs);
void read_sequence_sets(InputOptions iopt, const char * rfile, const char * ifile,
                         SequenceSet **rseq, SequenceSet **iseq);
double pdist(const char *a, const char *b, const int start, const int end,
             const int min_len);
int compute_distances(const SequenceSet *a, const char *seq,
                      const int start, const int end,
                      const int min_len,
                      double *pdistances);

int nucleotide2binary(const char *s, const int n, long unsigned int *b, long unsigned int *m, int *start, int *end);
SequenceSetB *read_aligned_sequencesB(const char *filename, const int len, const int num_seqs);
void read_sequence_setsB(InputOptions iopt, const char * rfile, const char * ifile,
                         SequenceSetB **rseq, SequenceSetB **iseq);
double pdistB(const long unsigned int *a, const long unsigned int *ma,
              const long unsigned int *b, const long unsigned int *mb,
              const int start, const int end, const int min_len);
int compute_distancesB(const SequenceSetB *a,
                       const long unsigned int *b, const long unsigned int *m,
                       const int start, const int end,
                       const int min_len,
                       double *pdistances);
distancesB_worker * init_distancesB_worker_pool(
    int nthread, const SequenceSetB *a, const int min_len, double *pdistances
);
int stop_distancesB_worker_pool(distancesB_worker* pool, int nthread);
int compute_distancesB_parallel(distancesB_worker * pool, int nthread,
                                const long unsigned int *b, const long unsigned int *m,
                                const int start, const int end);

/* routines_model.c */

Model *read_model(char *filename);
double **read_level_scalings(char *filename, int *num_levels);
int compute_cnode_probs_best2(const char *qid, TaxonomyNode *node, int nid, double prevprob, const Model *m, const double **scs, double pth, double rth, const double *pdistances);
int assign_with_discovery_best2(const char *qid, // name of the query
                                const int qi, // index of the query
                                Taxonomy *t, // full taxonomy
                                int nid, // id of node we are assigning within
                                const Model *m,
                                const double **scs,
                                const double *pdistances);
int print_model(Model *m);


#endif
