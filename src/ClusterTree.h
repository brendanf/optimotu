#ifndef OPTIMOTU_CLUSTERTREE_H_INCLUDED
#define OPTIMOTU_CLUSTERTREE_H_INCLUDED

#include "single_linkage.h"
#include "ClusterAlgorithm.h"
#include <deque>

#ifdef OPTIMOTU_R
#include <Rcpp.h>
#endif

// #define SINGLE_LINK_FULL_DEBUG
// #define SINGLE_LINK_DEBUG
// #define SINGLE_LINK_TEST


class ClusterTree : public SingleClusterAlgorithm {
protected:

  struct cluster {
    d_t min_d = NO_DIST;
    j_t id = NO_CLUST;
    uint32_t n_child = 0;
    cluster *parent = nullptr, *first_child = nullptr, *last_child = nullptr,
      *prev_sib = nullptr, *next_sib = nullptr;
    bool allocated = false;
    d_t max_d() {
      if (parent == nullptr) return NO_DIST;
      return parent->min_d;
    };
  };

   // tracks which clusters from the pool are currently used.
   std::deque<cluster*> freeclusters;
   std::vector<cluster> pool;
   cluster *const pool0, * const poolend, * const tip0, * const tipend, * const node0, * const nodeend;

   void delete_cluster(cluster * c) {
     freeclusters.push_back(c);
   };

   cluster * allocate_cluster() {
     cluster * c = freeclusters.front();
     freeclusters.pop_front();
     return c;
   };

   // assign all of the children of `csrc` to `cdest`
   // this DOES reassign the parent of all the children.
   void merge_children(cluster *cdest, cluster *csrc);

   // remove `child` from the linked child list belonging to `parent`
   // this DOES NOT reassign `c[child].parent`.
   void remove_child(cluster *parent, cluster *child);

   // add cluster `child` as a child of `parent`.
   // does NOT assign `child->parent`
   void add_child(cluster * parent, cluster * child);

   void shift_to_parent(cluster *& c, cluster *& cp) const;

#ifdef SINGLE_LINK_TEST
   void validate() const;

   std::string clust(const cluster * c) const {
     if (c) return std::to_string(c - my_pool);
     return "none";
   };

   std::string clust_id(const cluster * c) const {
     if (c) return std::to_string(c->id);
     return "none";
   };
#endif

#ifdef OPTIMOTU_R

   int hclust_ordering(cluster * top, int start, Rcpp::IntegerVector &order) const;

#endif

  void initialize() {
    j_t j;
    cluster *c;
    for (j = 0, c = tip0; c != tipend; j++, c++) {
      c->id = j; // id is smallest id of a descendent tip
      c->min_d = -1; // no minimum distance
      c->allocated = true;
    }
    for (; c < nodeend; c++) {
      freeclusters.push_back(c);
    }
  };

  void assign_ids();

  ClusterTree(SingleClusterAlgorithm * parent) :
    SingleClusterAlgorithm(parent),
    pool(2*n), pool0(pool.data()), poolend(pool0 + 2*n),
    tip0(pool0), tipend(tip0 + n),
    node0(tipend), nodeend(node0+n) {
    initialize();
  };
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

  ClusterTree(ClusterTree&& c) : SingleClusterAlgorithm(std::move(c)),
  freeclusters(std::move(c.freeclusters)),
  pool(std::move(c.pool)),
  pool0(pool.data()), poolend(pool0 + 2*n),
  tip0(pool0), tipend(tip0 + n),
  node0(tipend), nodeend(node0+n) {};

  ClusterTree * make_child() override;

  cluster* get_cluster(j_t j) const {
    if (j == NO_CLUST) return nullptr;
    return pool0 + j;
  };

   virtual void operator()(j_t seq1, j_t seq2, d_t i) override;

#ifdef OPTIMOTU_R
   void write_to_matrix(RcppParallel::RMatrix<int> &out) override;

   Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const override;
#endif

   // send consumer() pairwise distances to ensure it is up-to-date with this
   // clustering
   virtual void merge_into(DistanceConsumer &consumer) override;

   // send consumer() pairwise distances to ensure it is up-to-date with this
   // clustering
   virtual void merge_into(ClusterAlgorithm &consumer) override;

   // calculate the maximum distance between seq1 and seq2 which would actually
   // cause an update
   virtual double max_relevant(j_t seq1, j_t seq2) const override;
};

#endif //OPTIMOTU_CLUSTERTREE_H_INCLUDED
