#ifndef OPTIMOTU_CLUSTERTREE_H_INCLUDED
#define OPTIMOTU_CLUSTERTREE_H_INCLUDED

#include "single_linkage.h"
#include "ClusterAlgorithm.h"
#include <deque>
#include <cstdint>

#ifdef OPTIMOTU_R
#include <Rcpp.h>
#endif // OPTIMOTU_R

// #define CLUSTER_TREE_VERBOSE_DEBUG
// #define CLUSTER_TREE_DEBUG
// #define CLUSTER_TREE_TEST
// #define CLUSTER_TREE_FULL_TEST


class ClusterTree : public SingleClusterAlgorithm {
protected:

  struct cluster {
    d_t min_d = NO_DIST;
    j_t id = NO_CLUST;
    std::uint32_t n_child = 0;
    cluster *parent = nullptr, *first_child = nullptr, *last_child = nullptr,
      *prev_sib = nullptr, *next_sib = nullptr;
    bool allocated = false;
    d_t max_d();
  };

#ifdef CLUSTER_TREE_TEST
  struct cluster_int {
    const d_t min_d;
    const j_t id;
    const std::uint32_t n_child;
    const std::int32_t self, parent, first_child, last_child, prev_sib, next_sib;
    const bool allocated;

    cluster_int(const cluster * const c, const cluster * const pool0);
  };

  friend std::ostream& operator<<(std::ostream &out, const ClusterTree::cluster_int &c);
#endif

   // tracks which clusters from the pool are currently used.
   std::deque<cluster*> freeclusters;
   std::vector<cluster> pool;
   cluster *const pool0, * const poolend, * const tip0, * const tipend, * const node0, * const nodeend;

#ifdef CLUSTER_TREE_TEST
   size_t step_count = 0;
   mutable std::vector<cluster_int> touched_clusters;
   j_t current_seq1, current_seq2;
   d_t current_i;
#endif

   // marks the cluster as not-allocated and adds it to freeclusters
   void delete_cluster(cluster * c);

   // marks a new cluster as allocated, removes if from freeclusters
   cluster * allocate_cluster();

   // assign all of the children of `csrc` to `cdest` and delete `csrc`
   // this DOES reassign the parent of all the children.
   void merge_clusters(cluster *cdest, cluster *csrc);

   // remove `child` from the linked child list belonging to `parent`
   // this DOES NOT reassign `child->parent` or delete `child`.
   void remove_child(cluster *parent, cluster *child);

   // add cluster `child` as a child of `parent`.
   // does NOT assign `child->parent`
   void add_child(cluster * parent, cluster * child);

   void shift_to_parent(cluster *& c, cluster *& cp) const;

#ifdef CLUSTER_TREE_TEST
   bool validate_cluster(cluster * c) const;

   void validate_all() const;

   void validate_touched() const;

   std::string clust(const cluster * c) const;

   std::string clust_id(const cluster * c) const;
#endif // CLUSTER_TREE_TEST

#ifdef OPTIMOTU_R

   int hclust_ordering(cluster * top, int start, Rcpp::IntegerVector &order) const;

#endif // OPTIMOTU_R

  void initialize();

  void assign_ids();

  ClusterTree(SingleClusterAlgorithm * parent);

public:
   ClusterTree(const DistanceConverter &dconv, const j_t n):
     SingleClusterAlgorithm(dconv, n),
     pool(2*n), pool0(pool.data()), poolend(pool0 + 2*n),
     tip0(pool0), tipend(tip0 + n),
     node0(tipend), nodeend(node0+n) {
     initialize();
   };

  ClusterTree(const DistanceConverter &dconv, init_matrix_t im):
    SingleClusterAlgorithm(dconv, im),
    pool(2*n), pool0(pool.data()), poolend(pool0 + 2*n),
    tip0(pool0), tipend(tip0 + n),
    node0(tipend), nodeend(node0+n) {
    initialize();
  };

  ClusterTree * make_child() override;

  cluster* get_cluster(j_t j) const {
    if (j == NO_CLUST) return nullptr;
    return pool0 + j;
  };

  virtual void operator()(j_t seq1, j_t seq2, d_t i, int thread = 0) override;

  void write_to_matrix(internal_matrix_t &out) override;

#ifdef OPTIMOTU_R
  Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;
#endif // OPTIMOTU_R

  // send consumer() pairwise distances to ensure it is up-to-date with this
  // clustering
  virtual void merge_into(DistanceConsumer &consumer) override;

// send consumer() pairwise distances to ensure it is up-to-date with this
// clustering
virtual void merge_into(ClusterAlgorithm &consumer) override;

// calculate the maximum distance between seq1 and seq2 which would actually
// cause an update
virtual double max_relevant(j_t seq1, j_t seq2, int thread = 0) const override;
};

#endif //OPTIMOTU_CLUSTERTREE_H_INCLUDED
