#include <Rcpp.h>
#include <fstream>
#include <cmath>
#include <algorithm>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
#include "single_linkage_pool.hpp"

// #define SINGLE_LINK_FULL_DEBUG
// #define SINGLE_LINK_DEBUG
// #define SINGLE_LINK_TEST

// type of "j", which is the index of a cluster
typedef std::uint32_t j_t;
// type of "d", which is the index corresponding to a distance
typedef std::int32_t d_t;

struct cluster {
   d_t min_d = NO_DIST;
   j_t id = NO_CLUST, parent = NO_CLUST, n_child = 0,
      first_child = NO_CLUST, last_child = NO_CLUST,
      prev_sib = NO_CLUST, next_sib = NO_CLUST;
   bool allocated = false;
   d_t max_d(cluster *c) {
      if (parent == NO_CLUST) return NO_DIST;
      return c[parent].min_d;
   };
};

// tracks which clusters from the pool are currently used.
// note that it never actually interacts with the clusters.
struct cluster_pool {
   j_t *freeclusters;
   const int m;
   j_t n, first_free, last_free, n_free;
   cluster *my_pool;
   cluster_pool(int m, j_t n): m(m), n(n), n_free(n) {
      my_pool = new cluster[2*n];
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << "initializing first " << n << " of " << (2*n) <<
         " clusters" << std::endl;
#endif
      j_t j;
      cluster *c;
      for (j = 0, c = my_pool; j < n; j++, c++) {
         c->id = j; // id is smallest id of a descendent tip
         c->min_d = -1; // no minimum distance
         c->allocated = true;
      }
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << "finished cluster initialization" << std::endl;
#endif
      freeclusters = new j_t[n];
      for (j = 0; j < n; j++) {
         freeclusters[j] = n + j;
      }
      first_free = 0;
      last_free = n - 1;
   };

   ~cluster_pool() {
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << " in cluster_pool destructor" << std::endl;
      Rcpp::Rcout << "  freeing freeclusters at " << freeclusters << "..." << std::endl;
#endif
      delete[] freeclusters;
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << "  freeing my_pool at " << my_pool << "..." << std::endl;
#endif
      delete[] my_pool;
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << " leaving cluster_pool destructor" << std::endl;
#endif
   };

   cluster* get_cluster(j_t j) {
      if (j == NO_CLUST) return nullptr;
      return my_pool + j;
   };

   d_t max_d(j_t j) {
      if (j == NO_CLUST) return NO_DIST;
      return my_pool[j].max_d(my_pool);
   };

   d_t max_d(cluster * c) {
      if (c == nullptr) return NO_DIST;
      return c->max_d(my_pool);
   };

   void delete_cluster(j_t j) {
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << "-deallocating cluster " << j << "..." << std::endl;
#endif
      last_free++;
      if (last_free >= n) last_free = 0;
      freeclusters[last_free] = j;
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << " -last free cluster is " <<
         freeclusters[last_free] << " at position " <<
         last_free << "\n - " << n_free << " free clusters remain" << std::endl;
#endif
   };

   j_t allocate_cluster() {
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << "-allocating new cluster..." << std::endl;
#endif
      j_t j = freeclusters[first_free];
      first_free++;
      if (first_free >= n) first_free = 0;
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << " -allocated new cluster " << j <<
         "\n -first free cluster is " << freeclusters[first_free] <<
         " at position " << first_free << "\n - " << n_free <<
         " free clusters remain" << std::endl;
#endif
      return j;
   };
   // assign all of the children of `c[src]` to `c[dest]`
   // this DOES reassign the parent of all the children.
   void merge_children(j_t dest, j_t src) {
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << "-merging children of cluster " << src <<
         " into cluster " << dest <<
            " (clust array at " << my_pool << ")" << std::endl;
#endif
      cluster *csrc, *cdest, *c1, *c2;
      // CASE 1: the source is null
      // nothing to do
      if (src == NO_CLUST) {
#ifdef SINGLE_LINK_DEBUG
         Rcpp::Rcout << " -finished merging children of cluster " << src <<
            " into cluster " << dest << " (no-op)" << std::endl;
#endif
         return;
      }
      csrc = my_pool + src;
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << " - source cluster " << src << ": n_child=" << csrc->n_child <<
         ", first_child=" << csrc->first_child << ", last_child=" <<
            csrc->last_child << std::endl;
#endif
      // Reassign the parents
      // this is the slow part O(n)
      auto src_child = csrc->first_child;
      if (dest == NO_CLUST) {
#ifdef SINGLE_LINK_DEBUG
         Rcpp::Rcout << " - destination cluster " << dest <<
            " is null\n - removing children from source cluster..." << std::endl;
#endif
         while (src_child != NO_CLUST) {
#ifdef SINGLE_LINK_FULL_DEBUG
            Rcpp::Rcout << "  - removing child " << src_child << std::endl;
#endif
            c1 = my_pool + src_child;
            c1->parent = NO_CLUST;
            c1->prev_sib = NO_CLUST;
            c1->next_sib = NO_CLUST;
            src_child = c1->next_sib;
            csrc->n_child--;
         }
         csrc->first_child = NO_CLUST;
         csrc->last_child = NO_CLUST;
#ifdef SINGLE_LINK_DEBUG
         Rcpp::Rcout <<
            "  -finished removing children.\n  - source cluster " << src <<
               ": n_child=" << csrc->n_child <<
                  ", first_child=" << csrc->first_child << ", last_child=" <<
                     csrc->last_child << std::endl;
#endif
      } else {
         cdest = my_pool + dest;
#ifdef SINGLE_LINK_DEBUG
         Rcpp::Rcout << "  - destination cluster " << dest << ": n_child=" << cdest->n_child <<
            ", first_child=" << cdest->first_child << ", last_child=" <<
               cdest->last_child << "\n - merging..." << std::endl;
#endif
         // both clusters exist, so do the splice
         auto dest_child = cdest->last_child;
         if (dest_child != NO_CLUST && src_child != NO_CLUST) {
            c1 = my_pool + src_child;
            c2 = my_pool + dest_child;
            c2->next_sib = src_child;
            c1->prev_sib = dest_child;
         } else {
            cdest->first_child = src_child;
         }
         if (csrc->last_child != NO_CLUST) {
            cdest->last_child = csrc->last_child;
         }
         // reassign parents
         while (src_child != NO_CLUST) {
#ifdef SINGLE_LINK_FULL_DEBUG
            Rcpp::Rcout << "  - transferring child " << src_child << std::endl;
#endif
            c1 = my_pool + src_child;
            c1->parent = dest;
            csrc->n_child--;
            cdest->n_child++;
            src_child = c1->next_sib;
         }
         csrc->first_child = NO_CLUST;
         csrc->last_child = NO_CLUST;
#ifdef SINGLE_LINK_DEBUG
         Rcpp::Rcout <<
            "  -finished transferring children.\n  - source cluster " << src <<
               ": n_child=" << csrc->n_child <<
                  ", first_child=" << csrc->first_child << ", last_child=" <<
                     csrc->last_child << "\n  - destination cluster " << dest <<
                        ": n_child=" << cdest->n_child <<
                           ", first_child=" << cdest->first_child << ", last_child=" <<
                              cdest->last_child <<std::endl;
#endif
      }
      return;
   };

   // remove `child` from the linked child list belonging to `parent`
   // this DOES NOT reassign `c[child].parent`.
   void remove_child(j_t parent, j_t child) {
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << "-removing child " << child << " from cluster " << parent <<
         " (clust array at " << my_pool << ")" << std::endl;
#endif
      if (parent == NO_CLUST) {
#ifdef SINGLE_LINK_FULL_DEBUG
         Rcpp::Rcout << " -finished removing child " << child << " from cluster " <<
            parent << " (no-op)" << std::endl;
#endif
         return;
      }
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << "  -parent " << parent << ": n_child=" << my_pool[parent].n_child <<
         ", first_child=" << my_pool[parent].first_child << ", last_child=" << my_pool[parent].last_child <<
            "\n  -child " << child << ": prev_sib=" << my_pool[child].prev_sib <<
               ", next_sib=" << my_pool[child].next_sib << std::endl;
#endif

      auto prev = my_pool[child].prev_sib;
      auto next = my_pool[child].next_sib;
      if (prev == NO_CLUST) {
         my_pool[parent].first_child = next;
         if (next == NO_CLUST) {
            my_pool[parent].last_child = NO_CLUST;
            return;
         } else {
            my_pool[next].prev_sib = NO_CLUST;
            my_pool[child].next_sib = NO_CLUST;
         }
      } else {
         my_pool[prev].next_sib = next;
         if (next == NO_CLUST) {
            my_pool[parent].last_child = prev;
         } else {
            my_pool[next].prev_sib = prev;
            my_pool[child].next_sib = NO_CLUST;
         }
         my_pool[child].prev_sib = NO_CLUST;
      }
      my_pool[parent].n_child--;

#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << "  -finished removing child\n  -parent " << parent <<
         ": n_child=" << my_pool[parent].n_child <<
            ", first_child=" << my_pool[parent].first_child << ", last_child=" << my_pool[parent].last_child <<
               "\n  -child " << child << ": prev_sib=" << my_pool[child].prev_sib <<
                  ", next_sib=" << my_pool[child].next_sib << std::endl;
#endif
      return;
   };

   // add cluster `c[child]` as a child of `c[parent]`.
   // does NOT assign `c[child].parent`
   void add_child(j_t parent, j_t child) {
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << "-adding child " << child << " to cluster " << parent <<
         " (clust array at " << my_pool << ")" << std::endl;
#endif
      if (parent == NO_CLUST) {
#ifdef SINGLE_LINK_FULL_DEBUG
         Rcpp::Rcout << " -finished adding child " << child << " to cluster " <<
            parent << " (no-op)" << std::endl;
#endif
         return;
      }
      cluster *c_parent = my_pool + parent;
      cluster *c_child = my_pool + child;
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << "  -parent " << parent << ": n_child=" << c_parent->n_child <<
         ", first_child=" << c_parent->first_child << ", last_child=" << c_parent->last_child <<
            "\n  -child " << child << ": prev_sib=" << c_child->prev_sib <<
               ", next_sib=" << c_child->next_sib << std::endl;
#endif
      if (c_parent->last_child == NO_CLUST) {
         c_parent->first_child = c_parent->last_child = child;
         c_parent->n_child++;
#ifdef SINGLE_LINK_FULL_DEBUG
         Rcpp::Rcout << " -finished adding only child " << child << " to cluster " <<
            parent << "\n  -parent " << parent << ": n_child=" << c_parent->n_child <<
               ", first_child=" << c_parent->first_child << ", last_child=" << c_parent->last_child <<
                  "\n  -child " << child << ": prev_sib=" << my_pool[child].prev_sib <<
                     ", next_sib=" << c_child->next_sib << std::endl;
#endif
         return;
      }

      c_child->prev_sib = c_parent->last_child;
      my_pool[c_parent->last_child].next_sib = child;
      c_child->next_sib = NO_CLUST;
      c_parent->last_child = child;
      c_parent->n_child++;
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << " -finished adding child " << child << " to cluster " <<
         parent << "\n  -parent " << parent << ": n_child=" << c_parent->n_child <<
            ", first_child=" << c_parent->first_child << ", last_child=" << c_parent->last_child <<
               "\n  -child " << child << ": prev_sib=" << c_child->prev_sib <<
                  ", next_sib=" << c_child->next_sib << std::endl;
#endif
      return;
   }

   void shift_to_parent(j_t &j, cluster *& c, j_t &jp, cluster *& cp) {
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << "-shifting from child " << j << " to parent " << jp <<
         " (clust array at " << my_pool << ")" << std::endl;
#endif
      j = jp;
      c = cp;
      if (jp == NO_CLUST) {
#ifdef SINGLE_LINK_FULL_DEBUG
         Rcpp::Rcout << " -finished shifting to parent " << j << " (NO_CLUST)" <<
            std::endl;
#endif
         return;
      }
      jp = cp->parent;
      if (jp != NO_CLUST) cp = my_pool + jp; else cp = nullptr;
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << " -finished shifting to parent " << j << std::endl;
#endif
      return;
   }

#ifdef SINGLE_LINK_TEST
   void validate() {
      bool err = false;
      j_t j, jp;
      cluster *c, *cp;
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcerr << "validating..." << std::endl;
#endif
      for (j_t i = 0; i < n; i++) {
         j = i;
         c = my_pool + j;
         jp = c->parent;
         if (jp != NO_CLUST) {
            cp = my_pool + jp;
         } else {
            cp = nullptr;
         }
         while (true) {
            // Rcpp::Rcout << "cluster " << j <<
            //    "\n id: " << c->id <<
            //    "\n min_d: " << c->min_d <<
            //    "\n max_d: " << c->max_d(clust) <<
            //    "\n parent: " << c->parent <<
            //    "\n first_child: " << c->first_child <<
            //    "\n next_sib: " << c->next_sib <<
            //    "\n n_child: " << c->n_child << std::endl;
            if (c->min_d < -1 || (c->min_d >= m && c->min_d != NO_DIST)) {
               Rcpp::Rcerr << "validation error: cluster " << j <<
                  " has invalid min_d: " << c->min_d << std::endl;
               err = true;
            }
            std::int32_t max = c->max_d(my_pool);
            if (max < 0 || (max >= m && max != NO_DIST)) {
               Rcpp::Rcerr << "validation error: cluster " << j <<
                  " has invalid max_d: " << max << std::endl;
               err = true;
            }
            if (max < c->min_d) {
               Rcpp::Rcerr << "validation error: cluster " << j <<
                  " has invalid max_d: " << max <<
                     " which is less than its min_d: " << c->min_d << std::endl;
               err = true;
            }
            if (c->parent < n || (c->parent >= 2*n + 1 && c->parent != NO_CLUST)) {
               Rcpp::Rcerr << "validation error: cluster " << j <<
                  " has invalid parent: " << c->parent << std::endl;
               err = true;
            } else if (c->parent != NO_CLUST && !my_pool[c->parent].allocated) {
               Rcpp::Rcerr << "validation error: cluster " << j <<
                  "'s parent " << c->parent << " is not allocated " <<
                     std::endl;
               err = true;
            }
            if (c->n_child == 1) {
               Rcpp::Rcerr << "validation error: cluster " << j <<
                  " has exactly one child." << std::endl;
               err = true;
            }
            if (c->first_child < 0 || (c->first_child >= 2*n + 1 && c->first_child != NO_CLUST)) {
               Rcpp::Rcerr << "validation error: cluster " << j <<
                  " has invalid first_child: " << c->first_child << std::endl;
               err = true;
            }
            auto next = c->first_child;
            std::uint32_t kids = 0;
            while (next != NO_CLUST) {
               auto cnext = my_pool + next;
               kids++;
               if (cnext->parent != j) {
                  Rcpp::Rcerr << "validation error: cluster " << j << "'s child " <<
                     next << " instead has parent " << cnext->parent << std::endl;
                  err = true;
               }
               if (!cnext->allocated) {
                  Rcpp::Rcerr << "validation error: cluster " << j << "'s child " <<
                     next << " is not allocated " << std::endl;
                  err = true;
               }
               next = cnext->next_sib;
            }
            if (kids != c->n_child) {
               Rcpp::Rcerr << "validation error: cluster " << j <<
                  " has " << kids << " children but n_child=" << c->n_child << std::endl;
               err = true;
            }
            if (c->next_sib < 0 || (c->next_sib >= 2*n + 1 && c->next_sib != NO_CLUST)) {
               Rcpp::Rcerr << "validation error: cluster " << j <<
                  " has invalid next_sib: " << c->next_sib << std::endl;
               err = true;
            }
            if (c->next_sib != NO_CLUST && my_pool[c->next_sib].parent != c->parent) {
               Rcpp::Rcerr << "validation error: cluster " << j <<
                  " with parent " << c->parent << " has sibling " <<
                     c->next_sib << " which instead has parent " <<
                        my_pool[c->next_sib].parent << std::endl;
               err = true;
            }
            if (c->next_sib != NO_CLUST && !my_pool[c->next_sib].allocated) {
               Rcpp::Rcerr << "validation error: cluster " << j <<
                  " with parent " << c->parent << " has sibling " <<
                     c->next_sib << " which is not allocated " << std::endl;
               err = true;
            }
            if (c->next_sib != NO_CLUST && my_pool[c->next_sib].prev_sib != j) {
               Rcpp::Rcerr << "validation error: cluster " << j <<
                  " with parent " << c->parent << " has sibling " <<
                     c->next_sib << " whose previous sibling is " <<
                        my_pool[c->next_sib].prev_sib << std::endl;
               err = true;
            }
            if (c->next_sib == j ) {
               Rcpp::Rcerr << "validation error: cluster " << j <<
                  " is its own next sibling" << std::endl;
               err = true;
            }
            if (err) Rcpp::stop("found validation errors");
            if (max < m && c->prev_sib == NO_CLUST) shift_to_parent(j, c, jp, cp); else break;
         }
      }
   }
#endif

   void process(j_t seq1, j_t seq2, d_t i) {
      j_t j1 = seq1;
      cluster* c1 = get_cluster(seq1);
      j_t j1p = c1->parent;
      cluster* c1p = get_cluster(j1p);
      j_t j2 = seq2;
      cluster* c2 = get_cluster(seq2);
      j_t j2p = c2->parent;
      cluster* c2p = get_cluster(j2p);
      // search from the starting tips up the tree, until either the clusters
      // match (in which case nothing will happen) or both clusters are at the
      // correct height.
      d_t max_d1 = max_d(j1);
      d_t max_d2 = max_d(j2);
      while (j1 != j2 && (max_d1 <= i || max_d2 <= i)) {
         RcppThread::checkUserInterrupt();
         // shift to the parent of whichever cluster has the nearer parent
         if (max_d1 <= max_d2) {
            shift_to_parent(j1, c1, j1p, c1p);
            max_d1 = max_d(j1);
         } else {
            shift_to_parent(j2, c2, j2p, c2p);
            max_d2 = max_d(j2);
         }
      }
      bool changed = false;
      while (c1 && c2 && c1 != c2) {
         changed=true;
         /* we need to merge c1 and c2, and all their parents, starting at i and
          going up to m. */
#ifdef SINGLE_LINK_DEBUG
         Rcpp::Rcout << "clust: " << my_pool << "\n";
         Rcpp::Rcout << "j1:" << j1 << ", j2:" << j2 << ", i:" << i << "\n";
         Rcpp::Rcout << "min_d1:" << c1->min_d << ", min_d2:" << c2->min_d << "\n";
         Rcpp::Rcout << "max_d1:" << max_d(c1) << ", max_d2:" << max_d(c2) << std::endl;
#endif
         cluster *cnew;
         std::uint32_t jnew;
#ifdef SINGLE_LINK_DEBUG
         Rcpp::Rcout << "j1p:" << j1p << ", j2p:" << j2p << std::endl;
#endif
         if (i == c1->min_d) {
            // we don't need to create a new cluster for c1, we can just modify this
#ifdef SINGLE_LINK_DEBUG
            Rcpp::Rcout << "modifying first cluster: " << j1 << std::endl;
#endif
            jnew = j1;
            cnew = c1;
            if (i == c2->min_d) {
               remove_child(j2p, j2);
               // c2 is not needed anymore, so we need to assign its children to c1
               merge_children(j1, j2);
               // and mark it as unused
               c2->allocated = false;
               c2->parent = NO_CLUST;
               c2->first_child = NO_CLUST;
               c2->last_child = NO_CLUST;
               c2->next_sib = NO_CLUST;
               c2->prev_sib = NO_CLUST;
               delete_cluster(j2);
            } else {
               //add c2 as a child of cnew (aka c1)
               remove_child(j2p, j2);
               add_child(jnew, j2);
               c2->parent = jnew;
            }
         } else if (i == c2->min_d) {
            // we don't need to define a new cluster, we can re-use c2
#ifdef SINGLE_LINK_DEBUG
            Rcpp::Rcout << "modifying second cluster: " << j2 << std::endl;
#endif
            jnew = j2;
            cnew = c2;
            remove_child(j1p, j1);
            add_child(jnew, j1);
            c1->parent = jnew;
         } else if (c1p && c1p == c2p && c1p->n_child == 2) {
#ifdef SINGLE_LINK_DEBUG
            Rcpp::Rcout << "modifying common parent cluster: " << j1p << std::endl;
#endif
            // the two clusters to be merged have the same parent, and they are
            // the only children of that parent.
            // thus we can move the minimum distance of the parent rather than
            // creating a new cluster.
            cnew = c1p;
            cnew->min_d = i;
            shift_to_parent(j1, c1, j1p, c1p);
            shift_to_parent(j2, c2, j2p, c2p);
         } else {
            //create a new cluster as parent of current c1 and c2
            jnew = allocate_cluster();
            cnew = get_cluster(jnew);
            cnew->allocated = true;
            cnew->min_d = i;
            cnew->first_child = NO_CLUST;
            cnew->next_sib = NO_CLUST;
            cnew->parent = NO_CLUST;
            cnew->n_child = 0;
            remove_child(j1p, j1);
            add_child(jnew, j1); //sets cnew->n_child
            c1->parent = jnew;
            remove_child(j2p, j2);
            add_child(jnew, j2); //sets cnew->n_child
            c2->parent = jnew;
#ifdef SINGLE_LINK_DEBUG
            Rcpp::Rcout << "finished initializing cluster " << jnew << std::endl;
#endif
         }
         if (c1p && (c2p == nullptr || c2p->min_d >= c1p->min_d)) {
#ifdef SINGLE_LINK_DEBUG
            Rcpp::Rcout << "parent of new/modified cluster " << jnew <<
               " should be " << j1p << std::endl;
#endif
            if (cnew->parent != j1p) {
               remove_child(cnew->parent, jnew);
               cnew->parent = j1p;
               add_child(j1p, jnew);
            }
         } else {
            // either:
            // a) both are null
            // b) only c2p is non-null
            // c) both are non-null but c1p > c2p
#ifdef SINGLE_LINK_DEBUG
            Rcpp::Rcout << "parent of new/modified cluster " << jnew <<
               " should be " << j2p << std::endl;
#endif
            if (cnew->parent != j2p) {
               remove_child(cnew->parent, jnew);
               cnew->parent = j2p;
               add_child(j2p, jnew);
            }
         }
         // these will shift to j1p and c1p *even if the actual parent has changed*
         shift_to_parent(j1, c1, j1p, c1p);
         shift_to_parent(j2, c2, j2p, c2p);
         std::int32_t min_d1, min_d2;
         if (c1) {
            max_d1 = max_d(c1);
            min_d1 = c1->min_d;
         } else {
            max_d1 = min_d1 = NO_DIST;
         }
         if (c2) {
            max_d2 = max_d(c2);
            min_d2 = c2->min_d;
         } else {
            max_d2 = min_d2 = NO_DIST;
         }

         while (j1 != j2 && (c1p || c2p) && (max_d1 <= min_d2 || max_d2 <= min_d1)) {
            RcppThread::checkUserInterrupt();
#ifdef SINGLE_LINK_DEBUG
            // Rcpp::Rcout << "climbing tree: min_d1=" << min_d1 <<
            //    ", max_d1=" << max_d1 << ", min_d2=" << min_d2 <<
            //       ", max_d2=" << max_d2 << std::endl;
#endif
            // shift to the parent of whichever cluster has the nearer parent
            if (c1p && max_d1 <= max_d2) {
               shift_to_parent(j1, c1, j1p, c1p);
               max_d1 = max_d(c1);
               min_d1 = c1->min_d;
            } else {
               shift_to_parent(j2, c2, j2p, c2p);
               max_d2 = max_d(c2);
               min_d2 = c2->min_d;
            }
         }
         if (min_d1 < min_d2) i = min_d2; else i = min_d1;
      }
      if (c1 && c1 == c2 && c1->n_child < 2) {
         remove_child(j1p, j1);
         merge_children(j1p, j1);
         c1->parent = NO_CLUST;
         c1->allocated = false;
         c1->first_child = NO_CLUST;
         c1->last_child = NO_CLUST;
         c1->next_sib = NO_CLUST;
         c1->prev_sib = NO_CLUST;
         delete_cluster(j1);
      }
#ifdef SINGLE_LINK_TEST
      validate();
#endif
   }
};

//' @export
// [[Rcpp::export]]
Rcpp::IntegerMatrix single_linkage_pool(
    const std::string file,
    const Rcpp::CharacterVector &seqnames,
    const float dmin,
    const float dmax,
    const float dstep
) {
  const j_t n = seqnames.size();
  const int m = (int) ceilf((dmax - dmin)/dstep) + 1;
  cluster *c, *c1, *c2, *c1p, *c2p;
  j_t j, j1, j2, j1p, j2p;
  d_t i;

  // keep track of which clusters are free
  // clusters 0 to n-1 are always used (they are the tips)
  // to start with n to 2n-2 are free.
  cluster_pool pool(m, n);

  std::ifstream infile(file);
  j_t seq1, seq2;
  float dist;
  while(infile >> seq1 >> seq2 >> dist) {
    if (seq1 == seq2) continue;
#ifdef SINGLE_LINK_DEBUG
    Rcpp::Rcout << "seq1: " << seq1 << ", seq2:" << seq2 << ", dist:" << dist << std::endl;
#endif
    i = std::max((int) ceilf((dist - dmin) / dstep), 0);
    if (i >= m) continue;
    pool.process(seq1, seq2, i);
  }

  Rcpp::IntegerMatrix out(m, n);
  std::size_t k = 0;
  for (std::uint32_t i = 0; i < n; i++) {
     j = i;
     c = pool.get_cluster(j);
     j1p = c->parent;
     c1p = pool.get_cluster(j1p);
     std::uint32_t i2 = 0;
     while (i2 < m) {
        std::uint32_t max = pool.max_d(c);
        if (m < max) max = m;
        if (c->id == NO_CLUST) c->id = i;
        while (i2 < max) {
           out[k++] = c->id;
           i2++;
        }
        if (i2 < m) pool.shift_to_parent(j, c, j1p, c1p);
     }
  }
  return out;
}

void process(cluster_pool *pool, j_t seq1, j_t seq2, d_t i) {
   // RcppThread::Rcout << "dummy process; pool at " << pool << ", seq1=" <<
   //    seq1 << ", seq2=" << seq2 << ", i=" << i << std::endl;
   pool->process(seq1, seq2, i);
}

void process_mod(cluster_pool **pool, std::map<j_t, j_t> *fwd_map,
                 j_t seq1, j_t seq2, d_t i, std::vector<j_t> &whichsets,
                 std::uint8_t shard, std::uint8_t base) {
   for (j_t pc : whichsets) {
      if (pc % base != shard) continue;
      pool[pc]->process(fwd_map[pc][seq1], fwd_map[pc][seq2], i);
   }
}

//' @export
// [[Rcpp::export]]
Rcpp::List single_linkage_multi(
      const std::string file,
      const Rcpp::CharacterVector &seqnames,
      const float dmin,
      const float dmax,
      const float dstep,
      const Rcpp::ListOf<Rcpp::CharacterVector> &preclust,
      const size_t threads=1
) {
   RcppThread::ThreadPool::globalInstance().setNumThreads(threads);
   const j_t n = seqnames.size();
   const int m = (int) ceilf((dmax - dmin)/dstep) + 1;
   std::vector<j_t> *precluster_key = new std::vector<j_t>[n];
   // map from universal index to index in the precluster
   std::map<j_t, j_t> *fwd_map = new std::map<j_t, j_t>[preclust.size()];
   // keep track of which clusters are free
   // clusters 0 to n-1 are always used (they are the tips)
   // to start with n to 2n-2 are free.
   cluster_pool **pool = new cluster_pool*[preclust.size()];
   for (int pc = 0; pc < preclust.size(); pc++) {
      Rcpp::IntegerVector which = Rcpp::match(preclust[pc], seqnames);
      which = Rcpp::sort_unique(which) - 1;
      for(int j : which) {
#ifdef SINGLE_LINK_FULL_DEBUG
         Rcpp::Rcout << "adding seq " << j << " (" << seqnames[j] <<
            ") to precluster " << pc <<
            " at position " << fwd_map[pc].size() << std::endl;
#endif
         precluster_key[j].push_back(pc);
         fwd_map[pc].insert({j, fwd_map[pc].size()});

#ifdef SINGLE_LINK_FULL_DEBUG
         Rcpp::Rcout << "seq " << j << " now found in " <<
            precluster_key[j].size() << " preclusters \n precluster " << pc <<
               " now has " << fwd_map[pc].size() << "sequences" << std::endl;
#endif
      }
      pool[pc] = new cluster_pool(m, fwd_map[pc].size());
   }
   cluster *c, *c1, *c2, *c1p, *c2p;
   j_t j, j1, j2, j1p, j2p;
   d_t i;

   std::ifstream infile(file);
   j_t seq1, seq2;
   float dist;
   bool did_parallel = false;
   while(infile >> seq1 >> seq2 >> dist) {
      if (seq1 == seq2) continue;
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << "seq1: " << seq1 << ", seq2:" << seq2 << ", dist:" << dist << std::endl;
#endif
      i = std::max((int) ceilf((dist - dmin) / dstep), 0);
      if (i >= m) continue;
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << "seq " << seq1 << " present in preclusters:";
      for (auto pc : precluster_key[seq1]) Rcpp::Rcout << " " << pc;
      Rcpp::Rcout << std::endl;
      Rcpp::Rcout << "seq " << seq2 << " present in preclusters:";
      for (auto pc : precluster_key[seq2]) Rcpp::Rcout << " " << pc;
      Rcpp::Rcout << std::endl;
#endif

      std::vector<j_t> whichsets;
      whichsets.reserve(std::max(precluster_key[seq1].size(),
                                 precluster_key[seq2].size()));
      std::set_intersection(
         precluster_key[seq1].begin(),
         precluster_key[seq1].end(),
         precluster_key[seq2].begin(),
         precluster_key[seq2].end(),
         std::back_inserter(whichsets)
      );
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << "both seqs present in preclusters:";
      for (auto pc : whichsets) Rcpp::Rcout << " " << pc;
      Rcpp::Rcout << std::endl;
#endif
      if (did_parallel) {
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << "waiting for thread pool..." << std::endl;
#endif
      RcppThread::ThreadPool::globalInstance().wait();
      did_parallel = false;
      }
      if (threads > 1 && whichsets.size() > 1) {
         // RcppThread::parallelForEach(whichsets, [&pool, &fwd_map, seq1, seq2, i] (j_t pc) {
         //    pool[pc]->process(fwd_map[pc][seq1], fwd_map[pc][seq2], i);
         // });
         did_parallel = true;
         for (std::uint8_t k = 0; k < threads; k++) {
            RcppThread::ThreadPool::globalInstance().push(
                  process_mod,
               pool, fwd_map, seq1, seq2, i, whichsets, k, threads
            );
         }
//          for (j_t pc : whichsets) {
//             cluster_pool *pc_pool = pool[pc];
//             j_t pc_seq1 = fwd_map[pc][seq1];
//             j_t pc_seq2 = fwd_map[pc][seq2];
//             if (pc_pool->n >= 1000) {
//                did_parallel = true;
// #ifdef SINGLE_LINK_DEBUG
//                Rcpp::Rcout << "pushing thread to process seqs " << pc_seq1 << " and " <<
//                   pc_seq2 << " at threshold " << i << " for pool " << pc << std::endl;
// #endif
//                RcppThread::ThreadPool::globalInstance().push(
//                      process,
//                      pc_pool,
//                      pc_seq1,
//                      pc_seq2,
//                      i
//                );
//             } else {
//                process(pc_pool, pc_seq1, pc_seq2, i);
//             }
//          }
      } else {
         for (j_t pc : whichsets) {
            pool[pc]->process(fwd_map[pc][seq1], fwd_map[pc][seq2], i);
         }
      }
   }
#ifdef SINGLE_LINK_DEBUG
   Rcpp::Rcout << "waiting for thread pool..." << std::endl;
#endif
   RcppThread::ThreadPool::globalInstance().wait();

   Rcpp::List out(preclust.size());
   for(j_t pc = 0; pc < preclust.size(); pc++) {
      size_t pc_len = fwd_map[pc].size();
      Rcpp::IntegerMatrix outm(m, pc_len);

      std::size_t k = 0;
      for (j_t i = 0; i < pc_len; i++) {
         j = i;
         c = pool[pc]->get_cluster(j);
         j1p = c->parent;
         c1p = pool[pc]->get_cluster(j1p);
         j_t i2 = 0;
         while (i2 < m) {
            d_t max = pool[pc]->max_d(c);
            if (m < max) max = m;
            if (c->id == NO_CLUST) c->id = i;
            while (i2 < max) {
               outm[k++] = c->id;
               i2++;
            }
            if (i2 < m) pool[pc]->shift_to_parent(j, c, j1p, c1p);
         }
      }
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << " freeing pool " << pc << " at " << pool[pc] << std::endl;
#endif
      delete pool[pc];
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << " creating character names vector of size " << pc_len << std::endl;
#endif
      Rcpp::CharacterVector cn(pc_len);
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << " filling character names vector..." << std::endl;
#endif
      for (auto p : fwd_map[pc]) {
#ifdef SINGLE_LINK_FULL_DEBUG
         Rcpp::Rcout << " assigning seqname " << p.first <<
            " (" << seqnames[p.first] <<
            ") to column name " << p.second << std::endl;
#endif
         cn[p.second] = seqnames[p.first];
      }
      Rcpp::colnames(outm) = cn;
      out[pc] = outm;
   }
#ifdef SINGLE_LINK_DEBUG
   Rcpp::Rcout << " freeing pools..." << std::endl;
#endif
   delete[] pool;
#ifdef SINGLE_LINK_DEBUG
   Rcpp::Rcout << " freeing keys..." << std::endl;
#endif
   delete[] precluster_key;
#ifdef SINGLE_LINK_DEBUG
   Rcpp::Rcout << " freeing maps..." << std::endl;
#endif
   delete[] fwd_map;
   return out;
}
