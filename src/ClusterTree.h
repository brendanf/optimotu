#ifndef _CLUSTER_TREE_
#define _CLUSTER_TREE_
#include "single_linkage.h"
#include "ClusterAlgorithm.h"
#include <deque>

#ifdef OPTIMOTU_R
#include <Rcpp.h>
#endif

// #define SINGLE_LINK_FULL_DEBUG
// #define SINGLE_LINK_DEBUG
// #define SINGLE_LINK_TEST

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


class ClusterTree : public ClusterAlgorithm {
protected:
   // tracks which clusters from the pool are currently used.
   std::deque<cluster*> freeclusters;
   cluster *my_pool;

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
#endif

#ifdef OPTIMOTU_R

   int hclust_ordering(cluster * top, int start, Rcpp::IntegerVector &order) const;

#endif

  void initialize() {
    my_pool = new cluster[2*n];
    j_t j;
    cluster *c;
    for (j = 0, c = my_pool; j < n; j++, c++) {
      c->id = j; // id is smallest id of a descendent tip
      c->min_d = -1; // no minimum distance
      c->allocated = true;
    }
    for (; c < my_pool + 2*n; c++) {
      freeclusters.push_back(c);
    }
  };

  ClusterTree(ClusterAlgorithm * parent) :
    ClusterAlgorithm(parent) {
    initialize();
  };
public:
   ClusterTree(const DistanceConverter &dconv, const j_t n, const d_t m):
     ClusterAlgorithm(dconv, n, m) {
     initialize();
   };

  ClusterTree(ClusterTree&& c) : ClusterAlgorithm(std::move(c)),
  freeclusters(std::move(c.freeclusters)),
  my_pool(c.my_pool) {
    c.my_pool = nullptr;
  };

  ClusterTree * make_child() override;

   ~ClusterTree() {
      delete[] my_pool;
   };

   cluster* get_cluster(j_t j) const {
      if (j == NO_CLUST) return nullptr;
      return my_pool + j;
   };

   virtual void operator()(j_t seq1, j_t seq2, d_t i) override;

#ifdef OPTIMOTU_R
   void write_to_matrix(RcppParallel::RMatrix<int> &out) override;

   Rcpp::List as_hclust(const Rcpp::CharacterVector &seqnames) const;
#endif

   // send consumer() pairwise distances to ensure it is up-to-date with this
   // clustering
   virtual void merge_into(DistanceConsumer &consumer) const override;

   // send consumer() pairwise distances to ensure it is up-to-date with this
   // clustering
   virtual void merge_into(ClusterAlgorithm &consumer) const override;

   // calculate the maximum distance between seq1 and seq2 which would actually
   // cause an update
   virtual double max_relevant(j_t seq1, j_t seq2) const override;
};

#endif
