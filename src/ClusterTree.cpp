#include "ClusterTree.h"
#include "MultipleClusterAlgorithm.h"

#include <fstream>
#include <cmath>
#include <algorithm>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

void ClusterTree::operator()(j_t seq1, j_t seq2, d_t i, int thread) {
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  cluster* c1 = get_cluster(seq1);
  cluster* c1p = c1->parent;
  cluster* c2 = get_cluster(seq2);
  cluster* c2p = c2->parent;
  // search from the starting tips up the tree, until either the clusters
  // match (in which case nothing will happen) or both clusters are at the
  // correct height.
  d_t max_d1 = c1->max_d();
  d_t max_d2 = c2->max_d();
#ifdef CLUSTER_TREE_TEST
  ++step_count;
  current_seq1 = seq1;
  current_seq2 = seq2;
  current_i = i;
  touched_clusters.clear();
  touched_clusters.emplace_back(c1, this->pool0);
  if (c1p) touched_clusters.emplace_back(c1p, this->pool0);
  touched_clusters.emplace_back(c2, this->pool0);
  if (c2p) touched_clusters.emplace_back(c2p, this->pool0);
  int j = 0;
#endif
  while (c1 != c2 && (max_d1 <= i || max_d2 <= i)) {
    RcppThread::checkUserInterrupt();
    // shift to the parent of whichever cluster has the nearer parent
    if (max_d1 <= max_d2) {
      shift_to_parent(c1, c1p);
      max_d1 = c1->max_d();
    } else {
      shift_to_parent(c2, c2p);
      max_d2 = c2->max_d();
    }
#ifdef CLUSTER_TREE_TEST
    if (++j > 2*m) {
      OPTIMOTU_CERR << "infinite loop searching for clusters for sequences "
                    << seq1 << " and " << seq2
                    << " at threshold " << i << std::endl;
      validate_all();
      OPTIMOTU_STOP("infinite loop searching for clusters");
    }
#endif
  }
#ifdef CLUSTER_TREE_TEST
  j = 0;
#endif
  while (c1 && c2 && c1 != c2) {
    /* we need to merge c1 and c2, and all their parents, starting at i and
     going up to m. */
#ifdef CLUSTER_TREE_DEBUG
    Rcpp::Rcout << "j1:" << clust(c1)
                << ", j2:" << clust(c2)
                << ", i:" << i
                << std::endl;
    Rcpp::Rcout << "min_d1:" << c1->min_d
                << ", min_d2:" << c2->min_d
                << std::endl;
    Rcpp::Rcout << "max_d1:" << c1->max_d()
                << ", max_d2:" << c2->max_d()
                << std::endl;
#endif // CLUSTER_TREE_DEBUG
    cluster *cnew;
#ifdef CLUSTER_TREE_DEBUG
         Rcpp::Rcout << "j1p:" << clust(c1p)
                     << ", j2p:" << clust(c2p)
                     << std::endl;
#endif // CLUSTER_TREE_DEBUG
    if (i == c1->min_d) {
      // we don't need to create a new cluster for c1, we can just modify this
#ifdef CLUSTER_TREE_DEBUG
            Rcpp::Rcout << "modifying first cluster: " << clust(c1) << std::endl;
#endif // CLUSTER_TREE_DEBUG
      cnew = c1;
      if (i == c2->min_d) {
        remove_child(c2p, c2);
        // c2 is not needed anymore, so we need to assign its children to c1
        merge_children(c1, c2);
        // and mark it as unused
        c2->allocated = false;
        c2->parent = nullptr;
        c2->first_child = nullptr;
        c2->last_child = nullptr;
        c2->next_sib = nullptr;
        c2->prev_sib = nullptr;
        delete_cluster(c2);
      } else {
        //add c2 as a child of cnew (aka c1)
        remove_child(c2p, c2);
        add_child(cnew, c2);
        c2->parent = cnew;
      }
    } else if (i == c2->min_d) {
      // we don't need to define a new cluster, we can re-use c2
#ifdef CLUSTER_TREE_DEBUG
            Rcpp::Rcout << "modifying second cluster: " << clust(c2) << std::endl;
#endif // CLUSTER_TREE_DEBUG
      cnew = c2;
      remove_child(c1p, c1);
      add_child(cnew, c1);
      c1->parent = cnew;
    } else if (c1p && c1p == c2p && c1p->n_child == 2) {
#ifdef CLUSTER_TREE_DEBUG
            Rcpp::Rcout << "modifying common parent cluster: " << clust(c1p) << std::endl;
#endif // CLUSTER_TREE_DEBUG
      // the two clusters to be merged have the same parent, and they are
      // the only children of that parent.
      // thus we can move the minimum distance of the parent rather than
      // creating a new cluster.
      cnew = c1p;
      cnew->min_d = i;
      shift_to_parent(c1, c1p);
      shift_to_parent(c2, c2p);
    } else {
      //create a new cluster as parent of current c1 and c2
      cnew = allocate_cluster();
      cnew->allocated = true;
      cnew->min_d = i;
      cnew->first_child = nullptr;
      cnew->next_sib = nullptr;
      cnew->parent = nullptr;
      cnew->n_child = 0;
      remove_child(c1p, c1);
      add_child(cnew, c1); //sets cnew->n_child
      c1->parent = cnew;
      remove_child(c2p, c2);
      add_child(cnew, c2); //sets cnew->n_child
      c2->parent = cnew;
#ifdef CLUSTER_TREE_DEBUG
            Rcpp::Rcout << "finished initializing cluster " << clust(cnew) << std::endl;
#endif // CLUSTER_TREE_DEBUG
    }
    if (c1p && (c2p == nullptr || c2p->min_d >= c1p->min_d)) {
#ifdef CLUSTER_TREE_DEBUG
            Rcpp::Rcout << "parent of new/modified cluster " << clust(cnew)
                        << " should be " << clust(c1p)
                        << std::endl;
#endif // CLUSTER_TREE_DEBUG
      if (cnew->parent != c1p) {
        remove_child(cnew->parent, cnew);
        cnew->parent = c1p;
        add_child(c1p, cnew);
      }
    } else {
      // either:
      // a) both are null
      // b) only c2p is non-null
      // c) both are non-null but c1p > c2p
#ifdef CLUSTER_TREE_DEBUG
            Rcpp::Rcout << "parent of new/modified cluster " << clust(cnew)
                        << " should be " << clust(c2p)
                        << std::endl;
#endif // CLUSTER_TREE_DEBUG
      if (cnew->parent != c2p) {
        remove_child(cnew->parent, cnew);
        cnew->parent = c2p;
        add_child(c2p, cnew);
      }
    }
    // these will shift to j1p and c1p *even if the actual parent has changed*
    shift_to_parent(c1, c1p);
    shift_to_parent(c2, c2p);
    std::int32_t min_d1, min_d2;
    if (c1) {
      max_d1 = c1->max_d();
      min_d1 = c1->min_d;
    } else {
      max_d1 = min_d1 = NO_DIST;
    }
    if (c2) {
      max_d2 = c2->max_d();
      min_d2 = c2->min_d;
    } else {
      max_d2 = min_d2 = NO_DIST;
    }
#ifdef CLUSTER_TREE_TEST
    int k = 0;
#endif
    while (c1 != c2 && (c1p || c2p) && (max_d1 <= min_d2 || max_d2 <= min_d1)) {
      // shift to the parent of whichever cluster has the nearer parent
      if (c1p && max_d1 <= max_d2) {
        shift_to_parent(c1, c1p);
        max_d1 = c1->max_d();
        min_d1 = c1->min_d;
      } else {
        shift_to_parent(c2, c2p);
        max_d2 = c2->max_d();
        min_d2 = c2->min_d;
      }
#ifdef CLUSTER_TREE_TEST
      if (++k > m) {
        OPTIMOTU_CERR << "infinite loop recursing to parents "
                      << seq1 << " and " << seq2
                      << " at threshold" << i << std::endl;
        validate_touched();
        OPTIMOTU_STOP("infinite loop recursing to parents");
      }
#endif
    }
    if (min_d1 < min_d2) i = min_d2; else i = min_d1;
#ifdef CLUSTER_TREE_TEST
    if (++j > m) {
      OPTIMOTU_CERR << "infinite loop zipping subtrees for "
                    << seq1 << " and " << seq2
                    << " at threshold" << i << std::endl;
      validate_touched();
      OPTIMOTU_STOP("infinite loop zipping subtrees");
    }
#endif
  }
  if (c1 && c1 == c2 && c1->n_child < 2) {
    remove_child(c1p, c1);
    merge_children(c1p, c1);
    c1->parent = nullptr;
    c1->allocated = false;
    c1->first_child = nullptr;
    c1->last_child = nullptr;
    c1->next_sib = nullptr;
    c1->prev_sib = nullptr;
    delete_cluster(c1);
  }
#ifdef CLUSTER_TREE_TEST
  validate_touched();
#ifndef CLUSTER_TREE_FULL_TEST
  if (step_count % 100000 == 0) {
#endif //CLUSTER_TREE_FULL_TEST
  validate_all();
  OPTIMOTU_CERR << "finished " << step_count << " distances" << std::endl;
#ifndef CLUSTER_TREE_FULL_TEST
  }
#endif //CLUSTER_TREE_FULL_TEST
#endif // CLUSTER_TREE_TEST
}

d_t ClusterTree::cluster::max_d() {
  if (parent == nullptr) return NO_DIST;
  return parent->min_d;
}

#ifdef CLUSTER_TREE_TEST
ClusterTree::cluster_int::cluster_int(const cluster * const c, const cluster * const pool0):
  min_d(c->min_d), id(c->id), n_child(c->n_child),
  self(c - pool0),
  parent(c->parent ? c->parent - pool0 : -1),
  first_child(c->first_child ? c->first_child - pool0 : -1),
  last_child(c->last_child ? c->last_child - pool0 : -1),
  prev_sib(c->prev_sib ? c->prev_sib - pool0 : -1),
  next_sib(c->next_sib ? c->next_sib - pool0 : -1),
  allocated(c->allocated)
{}
#endif // CLUSTER_TREE_TEST

void ClusterTree::delete_cluster(cluster * c) {
  freeclusters.push_back(c);
}

ClusterTree::cluster * ClusterTree::allocate_cluster() {
  cluster * c = freeclusters.front();
  freeclusters.pop_front();
  return c;
}

void ClusterTree::initialize() {
  j_t j;
  cluster *c;
  for (j = 0, c = tip0; c != tipend; j++, c++) {
    c->id = j; // id is smallest id of a descendant tip
    c->min_d = -1; // no minimum distance
    c->allocated = true;
  }
  for (; c < nodeend; c++) {
    freeclusters.push_back(c);
  }
}

ClusterTree::ClusterTree(SingleClusterAlgorithm * parent) :
  SingleClusterAlgorithm(parent),
  pool(2*n), pool0(pool.data()), poolend(pool0 + 2*n),
  tip0(pool0), tipend(tip0 + n),
  node0(tipend), nodeend(node0+n) {
  initialize();
}

void ClusterTree::merge_children(cluster *cdest, cluster *csrc) {
#ifdef CLUSTER_TREE_DEBUG
   Rcpp::Rcout << "-merging children of cluster " << clust(csrc)
               << " into cluster " << clust(cdest)
               << std::endl;
#endif // CLUSTER_TREE_DEBUG
  // CASE 1: the source is null
  // nothing to do
  if (csrc == nullptr) {
#ifdef CLUSTER_TREE_DEBUG
         Rcpp::Rcout << " -finished merging children of cluster " << clust(csrc)
                     << " into cluster " << clust(cdest) << " (no-op)" << std::endl;
#endif // CLUSTER_TREE_DEBUG
    return;
  }
#ifdef CLUSTER_TREE_DEBUG
      Rcpp::Rcout << " - source cluster " << clust(csrc)
                  << ": n_child=" << csrc->n_child
                  << ", first_child=" << clust(csrc->first_child)
                  << ", last_child=" << clust(csrc->last_child)
                  << std::endl;
#endif // CLUSTER_TREE_DEBUG
  // Reassign the parents
  // this is the slow part O(n)
  auto src_child = csrc->first_child;
  if (cdest == nullptr) {
#ifdef CLUSTER_TREE_DEBUG
         Rcpp::Rcout << " - destination cluster " << clust(cdest)
                     << " is null" << std::endl
                     << " - removing children from source cluster..." << std::endl;
#endif // CLUSTER_TREE_DEBUG
    while (src_child) {
#ifdef CLUSTER_TREE_TEST
      touched_clusters.emplace_back(src_child, this->pool0);
      if (csrc->n_child == 0) {
        OPTIMOTU_CERR << "number of children larger than n_child in cluster "
                      << clust(csrc)
                      << " while merging with cluster "
                      << clust(cdest) << std::endl;
        validate_touched();
        OPTIMOTU_STOP("too many children");
      }
#endif // CLUSTER_TREE_TEST
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
            Rcpp::Rcout << "  - removing child " << clust(src_child) << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG
      src_child->parent = nullptr;
      src_child->prev_sib = nullptr;
      src_child->next_sib = nullptr;
      src_child = src_child->next_sib;
      csrc->n_child--;
    }
    csrc->first_child = nullptr;
    csrc->last_child = nullptr;
#ifdef CLUSTER_TREE_DEBUG
         Rcpp::Rcout << "  -finished removing children."
                     << std::endl
                     << "  - source cluster " << clust(csrc)
                     << ": n_child=" << csrc->n_child
                     << ", first_child=" << clust(csrc->first_child)
                     << ", last_child=" << clust(csrc->last_child)
                     << std::endl;
#endif // CLUSTER_TREE_DEBUG
  } else {
#ifdef CLUSTER_TREE_DEBUG
         Rcpp::Rcout << "  - destination cluster " << clust(cdest)
                     << ": n_child=" << cdest->n_child
                     << ", first_child=" << clust(cdest->first_child)
                     << ", last_child=" << clust(cdest->last_child)
                     << std::endl << " - merging..." << std::endl;
#endif // CLUSTER_TREE_DEBUG
    // both clusters exist, so do the splice
    auto dest_child = cdest->last_child;
    if (dest_child && src_child) {
#ifdef CLUSTER_TREE_TEST
      touched_clusters.emplace_back(dest_child, this->pool0);
#endif // CLUSTER_TREE_TEST
      dest_child->next_sib = src_child;
      src_child->prev_sib = dest_child;
    } else {
      cdest->first_child = src_child;
    }
    if (csrc->last_child) {
      cdest->last_child = csrc->last_child;
    }
    // reassign parents
    while (src_child) {
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
            Rcpp::Rcout << "  - transferring child " << clust(src_child) << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG
#ifdef CLUSTER_TREE_TEST
      touched_clusters.emplace_back(src_child, this->pool0);
      if (csrc->n_child == 0) {
        OPTIMOTU_CERR << "number of children larger than n_child in cluster "
                      << clust(csrc)
                      << " while merging with cluster "
                      << clust(cdest) << std::endl;
        validate_touched();
        OPTIMOTU_STOP("too many children");
      }
#endif // CLUSTER_TREE_TEST
      src_child->parent = cdest;
      csrc->n_child--;
      cdest->n_child++;
      src_child = src_child->next_sib;
    }
    csrc->first_child = nullptr;
    csrc->last_child = nullptr;
#ifdef CLUSTER_TREE_DEBUG
         Rcpp::Rcout << "  -finished transferring children."
                     << "  - source cluster " << clust(csrc)
                     << ": n_child=" << csrc->n_child
                     << ", first_child=" << clust(csrc->first_child)
                     << ", last_child=" << clust(csrc->last_child)
                     << std::endl
                     << "  - destination cluster " << clust(cdest)
                     << ": n_child=" << cdest->n_child
                     << ", first_child=" << clust(cdest->first_child)
                     << ", last_child=" << clust(cdest->last_child)
                     << std::endl;
#endif // CLUSTER_TREE_DEBUG
  }
  return;
}

void ClusterTree::remove_child(cluster *parent, cluster *child) {
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
      Rcpp::Rcout << "-removing child " << clust(child)
                  << " from cluster " << clust(parent)
                  << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG
  if (parent == nullptr) {
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
         Rcpp::Rcout << " -finished removing child " << clust(child)
                     << " from cluster " << clust(parent)
                     << " (no-op)" << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG
    return;
  }
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
      Rcpp::Rcout << "  -parent " << clust(parent)
                  << ": n_child=" << parent->n_child
                  << ", first_child=" << clust(parent->first_child)
                  << ", last_child=" << clust(parent->last_child)
                  << std::endl
                  << "  -child " << clust(child)
                  << ": prev_sib=" << clust(child->prev_sib)
                  << ", next_sib=" << clust(child->next_sib)
                  << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG

  auto prev = child->prev_sib;
  auto next = child->next_sib;
  if (prev == nullptr) {
    parent->first_child = next;
    if (next == nullptr) {
      parent->last_child = nullptr;
      return;
    } else {
      next->prev_sib = nullptr;
      child->next_sib = nullptr;
    }
  } else {
    prev->next_sib = next;
    if (next == nullptr) {
      parent->last_child = prev;
    } else {
      next->prev_sib = prev;
      child->next_sib = nullptr;
    }
    child->prev_sib = nullptr;
  }
  parent->n_child--;
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
      Rcpp::Rcout << "  -finished removing child" << std::endl
                  << "  -parent " << clust(parent)
                  << ": n_child=" << parent->n_child
                  << ", first_child=" << clust(parent->first_child)
                  << ", last_child=" << clust(parent->last_child)
                  << std::endl
                  << "  -child " << clust(child)
                  << ": prev_sib=" << clust(child->prev_sib)
                  << ", next_sib=" << clust(child->next_sib)
                  << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG
  return;
}

void ClusterTree::add_child(cluster * parent, cluster * child) {
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
      Rcpp::Rcout << "-adding child " << clust(child)
                  << " to cluster " << clust(parent)
                  << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG
  if (parent == nullptr) {
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
         Rcpp::Rcout << " -finished adding child " << clust(child)
                     << " to cluster " << clust(parent)
                     << " (no-op)" << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG
    return;
  }
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
      Rcpp::Rcout << "  -parent " << clust(parent)
                  << ": n_child=" << parent->n_child
                  << ", first_child=" << clust(parent->first_child)
                  << ", last_child=" << clust(parent->last_child)
                  << std::endl
                  << "  -child " << clust(child)
                  << ": prev_sib=" << clust(child->prev_sib)
                  << ", next_sib=" << clust(child->next_sib)
                  << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG
  if (parent->last_child == nullptr) {
#ifdef CLUSTER_TREE_TEST
    if (parent < this->node0) {
      OPTIMOTU_CERR << "attempting to add child " << clust(child)
                    << " to leaf node " << clust(parent)
                    << std::endl;
      OPTIMOTU_STOP("tried to add child to leaf node");
    }
#endif // CLUSTER_TREE_TEST
    parent->first_child = parent->last_child = child;
    parent->n_child++;
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
         Rcpp::Rcout << " -finished adding only child " << clust(child)
                     << " to cluster " << clust(parent)
                     << std::endl
                     << "  -parent " << clust(parent)
                     << ": n_child=" << parent->n_child
                     << ", first_child=" << clust(parent->first_child)
                     << ", last_child=" << clust(parent->last_child)
                     << std::endl
                     << "  -child " << clust(child)
                     << ": prev_sib=" << clust(child->prev_sib)
                     << ", next_sib=" << clust(child->next_sib)
                     << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG
    return;
  }

  child->prev_sib = parent->last_child;
  parent->last_child->next_sib = child;
  child->next_sib = nullptr;
  parent->last_child = child;
  parent->n_child++;
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
      Rcpp::Rcout << " -finished adding child " << clust(child)
                  << " to cluster " << clust(parent)
                  << std::endl
                  << "  -parent " << clust(parent)
                  << ": n_child=" << parent->n_child
                  << ", first_child=" << clust(parent->first_child)
                  << ", last_child=" << clust(parent->last_child)
                  << std::endl
                  << "  -child " << clust(child)
                  << ": prev_sib=" << clust(child->prev_sib)
                  << ", next_sib=" << clust(child->next_sib)
                  << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG
  return;
}

void ClusterTree::shift_to_parent(cluster *& c, cluster *& cp) const {
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
      Rcpp::Rcout << "-shifting from child " << clust(c)
                  << " to parent " << clust(cp)
                  << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG
  c = cp;
  if (cp == nullptr) {
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
         Rcpp::Rcout << " -finished shifting to parent " << clust(c)
                     << " (NO_CLUST)"
                     << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG
    return;
  }
  cp = cp->parent;
#ifdef CLUSTER_TREE_TEST
  if (cp != nullptr) {
    touched_clusters.emplace_back(cp, this->pool0);
  }
#endif // CLUSTER_TREE_TEST
#ifdef CLUSTER_TREE_VERBOSE_DEBUG
      Rcpp::Rcout << " -finished shifting to parent " << clust(c) << std::endl;
#endif // CLUSTER_TREE_VERBOSE_DEBUG
  return;
}

#ifdef OPTIMOTU_R
int ClusterTree::hclust_ordering(cluster * top, int start, Rcpp::IntegerVector &order) const {
  // Rcpp::Rcout << "ordering cluster " << top << " with start=" << start << std::endl;
  int i = start;
  cluster *c = top->first_child;
  if (c == nullptr) {
    // Rcpp::Rcout << i << ": " << top << std::endl;
    order[i++] = (top - this->pool0) + 1;
  } else {
    do {
      // Rcpp::Rcout << "getting cluster " << j << std::endl;
      if (c < this->tipend) {
        // Rcpp::Rcout << i << ": " << j << std::endl;
        order[i++] = c - this->pool0 + 1;
        // i++;
        // Rcpp::Rcout << "i=" << i << std::endl;
      } else {
        i = hclust_ordering(c, i, order);
      }
      c = c->next_sib;
    } while (c);
  }
  return i;
}
#endif // OPTIMOTU_R

void ClusterTree::assign_ids() {
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  for (j_t i = 0; i < n; ++i) {
    cluster * c = &pool[i];
    while(c->parent && c->parent->id > i) {
      c->parent->id = i;
      c = c->parent;
    }
  }
}

void ClusterTree::merge_into(DistanceConsumer &consumer) {
  this->assign_ids();
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  for (auto c = this->node0; c < this->nodeend; ++c) {
    if (c->allocated) {
      if (c->first_child) {
        cluster * next = c->first_child->next_sib;
        while (next) {
          consumer(c->first_child->id, next->id, dconv.inverse(c->min_d));
          next = next->next_sib;
        }
      }
    }
  }
}

void ClusterTree::merge_into(ClusterAlgorithm &consumer) {
  this->assign_ids();
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  for (auto c = this->node0; c < this->nodeend; ++c) {
    // OPTIMOTU_COUT << "merging cluster " << c - this->pool0
    //           << " (" << c
    //           << ") with ID " << c->id
    //           << std::endl;
    if (c->allocated) {
      // OPTIMOTU_COUT << " - cluster is allocated" << std::endl;
      if (c->first_child) {
        // OPTIMOTU_COUT << " - first child is cluster " << c->first_child - this->pool0
        //           << " (" << c->first_child
        //           << ") with ID " << c->first_child->id
        //           << std::endl;
        cluster * next = c->first_child->next_sib;
        while (next) {
          // OPTIMOTU_COUT << " - next child is cluster " << next - this->pool0
          //           << " (" << next
          //           << ") with ID " << next->id
          //           << std::endl;
          consumer(c->first_child->id, next->id, c->min_d);
          next = next->next_sib;
        }
      }
    }
  }
}

double ClusterTree::max_relevant(j_t seq1, j_t seq2, int thread) const {
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  cluster* c1 = this->get_cluster(seq1);
  cluster* c1p = c1->parent;
  cluster* c2 = this->get_cluster(seq2);
  cluster* c2p = c2->parent;
  // search from the starting tips up the tree, until either the clusters
  // match (in which case nothing will happen) or both clusters are at the
  // correct height.
  d_t max_d1 = c1->max_d();
  d_t max_d2 = c2->max_d();
  if (c1p == c2p && max_d1 == 0) return -1.0;
  while (c1 != c2) {
    RcppThread::checkUserInterrupt();
    // shift to the parent of whichever cluster has the nearer parent
    if (max_d1 <= max_d2 && c1p) {
      this->shift_to_parent(c1, c1p);
      max_d1 = c1->max_d();
    } else if (max_d2 <= max_d1 && c2p){
      this->shift_to_parent(c2, c2p);
      max_d2 = c2->max_d();
    } else {
      break;
    }
  }
  if (c1 == c2) return this->dconv.inverse(c1->min_d - 1);
  return this->dconv.inverse(this->m - 1);
}


void ClusterTree::write_to_matrix(internal_matrix_t &out) {
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  j_t j;
  cluster *c, *c1p;
  std::size_t k = 0;
  std::size_t n = out.ncol();
  std::size_t m = out.nrow();
  for (std::uint32_t i = 0; i < n; i++) {
    j = i;
    c = get_cluster(j);
    c1p = c->parent;
    j_t i2 = 0;
    while (i2 < m) {
      d_t max = c->max_d();
      if (m < (size_t)max) max = m;
      if (c->id == NO_CLUST) c->id = i;
      while (i2 < (size_t)max) {
        out[k++] = c->id;
        i2++;
      }
      if (i2 < m) shift_to_parent(c, c1p);
    }
  }
}

#ifdef OPTIMOTU_R
Rcpp::List ClusterTree::as_hclust(const Rcpp::CharacterVector &seqnames) const {
  std::shared_lock<std::shared_timed_mutex> lock(this->mutex);
  Rcpp::IntegerMatrix merge(this->n - 1, 2);
  Rcpp::NumericVector height(this->n-1);
  Rcpp::IntegerVector order(this->n);

  // Rcpp::Rcout << "making cluster/distance pairs" << std::endl;
  // depths and indices of the valid clusters
  std::vector<std::pair<d_t, cluster*>> cluster_dj;
  for (cluster *c = this->pool0; c < this->poolend; ++c) {
    // Rcpp::Rcout << "adding cluster " << c << " at depth " << c->min_d << std::endl;
    if (c->allocated) cluster_dj.emplace_back(c->min_d, c);
  }

  // Rcpp::Rcout << "sorting cluster/distance pairs" << std::endl;
  std::sort(cluster_dj.begin(), cluster_dj.end());
  // Rcpp::Rcout << "finished sorting cluster/distance pairs" << std::endl;
  size_t row = 0; // 1-indexed for R
  auto merge1 = merge.begin(); // first column of merge
  auto merge2 = merge1 + this->n - 1; // second column of merge
  auto h = height.begin();
  std::map<cluster*, int> row_key; // maps cluster to row number
  std::vector<cluster*> parentless; //keep track of clusters with no parent
  for (auto cli : cluster_dj) {
    // Rcpp::Rcout << "processing cluster " << cli.second << " at depth " << cli.first << std::endl;
    cluster *c = cli.second;
    cluster *child_c = c->first_child;
    cluster *next_c = child_c->next_sib;
    cluster *left = child_c;
    bool first_pass = TRUE;
    do {
      if (left - this->pool0 < this->n && first_pass) {
        // Rcpp::Rcout << "left child = " << -left << std::endl;
        *merge1 = this->pool0 - left - 1;
      } else if (first_pass) {
        // Rcpp::Rcout << "left child = " << row_key[child_j] << std::endl;
        *merge1 = row_key[child_c];
      } else {
        // Rcpp::Rcout << "left child = " << row << std::endl;
        *merge1 = row; //not incremented yet; this is the previous row
      }
      if (next_c - this->pool0 < this->n) {
        // Rcpp::Rcout << "right child = " << -(int)next_j << std::endl;
        *merge2 = (int)(this->pool0 - next_c - 1);
      } else {
        // Rcpp::Rcout << "right child = " << row_key[next_j] << std::endl;
        *merge2 = row_key[next_c];
      }
      // Rcpp::Rcout << "height = " << dconv.inverse(c->min_d) << std::endl;
      *h = dconv.inverse(c->min_d);
      // go to the next row
      merge2++;
      merge1++;
      h++;
      row++;
      first_pass = FALSE;
      child_c = next_c;
      next_c = child_c->next_sib;
    } while (next_c != nullptr);
    if (c->parent == nullptr) parentless.push_back(cli.second);
    // Rcpp::Rcout << "finished processing cluster " << cli.second << " at row " << row << std::endl;

    row_key[cli.second] = row;
  }
  // Rcpp::Rcout << "finished processing clusters" << std::endl;
  // find parentless singletons
  for (cluster * c = this->tip0; c < this->tipend; c++) {
    if (c->parent == nullptr) {
      parentless.push_back(c);
    }
  }
  if (parentless.size() > 1) {
    // Rcpp::Rcout << "adding dummy parent for parentless clusters" << std::endl;
    int left;
    if (parentless[0] < this->tipend) {
      left = (int)(this->pool0 - parentless[0] - 1);
    } else {
      left = row_key[parentless[0]];
    }
    for (size_t i = 1; i < parentless.size(); i++) {
      int right;
      if (parentless[i] < this->tipend) {
        right = (int)(this->pool0 - parentless[i] - 1);
      } else {
        right = row_key[parentless[i]];
      }
      // Rcpp::Rcout << "left child = " << left << std::endl;
      *(merge1++) = left;
      // Rcpp::Rcout << "right child = " << right << std::endl;
      *(merge2++) = right;
      *(h++) = 1;
      row++;
      left = row;
    }
    // Rcpp::Rcout << "finished adding parentless clusters" << std::endl;
  }
  j_t i = 0;
  for (cluster * c : parentless) {
    // Rcpp::Rcout << "calculating ordering for parentless cluster " << j << std::endl;
    i = hclust_ordering(c, i, order);
  }

  Rcpp::List out;
  out["merge"] = merge;
  out["height"] = height;
  out["order"] = order;
  out["labels"] = seqnames;
  out["method"] = "single";
  out.attr("class") = "hclust";
  return out;
}
#endif // OPTIMOTU_R

ClusterTree * ClusterTree::make_child() {
  std::unique_lock<std::shared_timed_mutex> lock(this->mutex);
  if (own_child) {
    auto child_ptr = new ClusterTree(this);
    auto child = std::unique_ptr<ClusterAlgorithm>(
      (ClusterAlgorithm*)child_ptr
    );
    this->children.push_back(std::move(child));
    return child_ptr;
  }
  this->own_child = true;
  return this;
}

#ifdef CLUSTER_TREE_TEST
std::string ClusterTree::clust(const cluster * c) const {
  if (c) return std::to_string(c - pool0);
  return "none";
}

std::string ClusterTree::clust_id(const cluster * c) const {
  if (c) return std::to_string(c->id);
  return "none";
}

inline std::ostream& operator<<(std::ostream &out, const ClusterTree::cluster_int &c)
{
  out << "  -" << std::endl
      << "    self: " << c.self << std::endl
      << "    allocated: " << c.allocated << std::endl
      << "    id: " << c.id << std::endl
      << "    min_d: " << c.min_d << std::endl
      << "    parent: " << c.parent << std::endl
      << "    first_child: " << c.first_child << std::endl
      << "    last_child: " << c.last_child << std::endl
      << "    n_child: " << c.n_child << std::endl
      << "    prev_sib: " << c.prev_sib << std::endl
      << "    next_sib: " << c.next_sib << std::endl;
  return out;
}

bool ClusterTree::validate_cluster(cluster * c) const {

  if (!c->allocated) return true;

  bool status = true;

  if (c->min_d < -1 || (c->min_d >= m && c->min_d != NO_DIST)) {
    OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                  << " has invalid min_d: " << c->min_d << std::endl;
    status = false;
  }
  std::int32_t max = c->max_d();
  if (max < 0 || (max >= m && max != NO_DIST)) {
    OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                  << " has invalid max_d: " << max << std::endl;
    status = false;
  }
  if (max < c->min_d) {
    OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                  << " has invalid max_d: " << max
                  << " which is less than its min_d: " << c->min_d << std::endl;
    status = false;
  }
  if (c->parent && (c->parent < this->node0 ||
      c->parent > this->nodeend)) {
    OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                  << " has invalid parent: " << clust(c->parent)
                  << std::endl;
    status = false;
  } else if (c->parent && !c->parent->allocated) {
    OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                  << "'s parent " << clust(c->parent) << " is not allocated "
                  << std::endl;
    status = false;
  }
  if (c->n_child == 1) {
    OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                  << " has exactly one child." << std::endl;
    status = false;
  }
  if (c -> first_child && (c->first_child < this->pool0
                             || c->first_child > this->poolend)) {
    OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                  << " has invalid first_child: " << clust(c->first_child) << std::endl;
    status = false;
  }
  auto next = c->first_child;
  std::uint32_t kids = 0;
  while (next) {
    kids++;
    if (next->parent != c) {
      OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                    << "'s child " << clust(next)
                    << " instead has parent " << clust(next->parent)
                    << std::endl;
      status = false;
    }
    if (!next->allocated) {
      OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                    << "'s child " << clust(next)
                    << " is not allocated " << std::endl;
      status = false;
    }
    next = next->next_sib;
  }
  if (kids != c->n_child) {
    OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                  << " has " << kids
                  << " children but n_child=" << c->n_child
                  << std::endl;
    status = false;
  }
  if (c -> next_sib && (c->next_sib < this->pool0
                          || c->next_sib > this->poolend)) {
    OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                  << " has invalid next_sib: " << clust(c->next_sib)
                  << std::endl;
    status = false;
  }
  if (c->next_sib && c->next_sib->parent != c->parent) {
    OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                  << " with parent " << clust(c->parent)
                  << " has sibling " <<  clust(c->next_sib)
                  << " which instead has parent " << clust(c->next_sib->parent)
                  << std::endl;
    status = false;
  }
  if (c->next_sib && !c->next_sib->allocated) {
    OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                  << " with parent " << clust(c->parent)
                  << " has sibling " << clust(c->next_sib)
                  << " which is not allocated " << std::endl;
    status = false;
  }
  if (c->next_sib && c->next_sib->prev_sib != c) {
    OPTIMOTU_CERR << "validation error: cluster " << c - this->pool0
                  << " with parent " << clust(c->parent)
                  << " has sibling " << clust(c->next_sib)
                  << " whose previous sibling is " << clust(c->next_sib->prev_sib)
                  << std::endl;
    status = false;
  }
  if (c->next_sib == c ) {
    OPTIMOTU_CERR << "validation error: cluster " << clust(c)
                  << " is its own next sibling" << std::endl;
    status = false;
  }
  return status;
}

void ClusterTree::validate_all() const {
  bool status = true;
#ifdef CLUSTER_TREE_DEBUG
  OPTIMOTU_CERR << "validating all clusters..." << std::endl;
#endif // CLUSTER_TREE_DEBUG
  for (cluster * c = this->pool0; c < this->poolend; ++c) {
    status &= validate_cluster(c);
  }
  if (!status) OPTIMOTU_STOP("found validation errors");
}

void ClusterTree::validate_touched() const {
  bool status = true;
#ifdef CLUSTER_TREE_DEBUG
  OPTIMOTU_CERR << "validating touched clusters..." << std::endl;
#endif //CLUSTER_TREE_DEBUG
  for (const cluster_int &c : touched_clusters) {
    status &= validate_cluster(this->pool0 + c.self);
  }
  if (!status) {
    OPTIMOTU_CERR << "error while processing pairwise distance "
                  << step_count
                  << " between seq1=" << this->current_seq1
                  << ", seq2=" << this->current_seq2
                  << ", i=" << this->current_i
                  << std::endl
                  << "---" << std::endl
                  << "clusters:" << std::endl;
    for (const cluster_int &c : touched_clusters) {
      OPTIMOTU_CERR << c;
    }
    OPTIMOTU_CERR << "---" << std::endl;
    OPTIMOTU_STOP("found validation errors");
  }

}


#endif // CLUSTER_TREE_TEST
