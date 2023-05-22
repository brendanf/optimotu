#include "ClusterTree.h"
#include "MultipleClusterAlgorithm.h"

#include <fstream>
#include <cmath>
#include <algorithm>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>

void ClusterTree::operator()(j_t seq1, j_t seq2, d_t i) {
  // std::lock_guard<std::mutex> lock(this->mutex);
  tbb::queuing_rw_mutex::scoped_lock lock(this->mutex);
  cluster* c1 = get_cluster(seq1);
  cluster* c1p = c1->parent;
  cluster* c2 = get_cluster(seq2);
  cluster* c2p = c2->parent;
  // search from the starting tips up the tree, until either the clusters
  // match (in which case nothing will happen) or both clusters are at the
  // correct height.
  d_t max_d1 = c1->max_d();
  d_t max_d2 = c2->max_d();
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
  }
  while (c1 && c2 && c1 != c2) {
    /* we need to merge c1 and c2, and all their parents, starting at i and
     going up to m. */
#ifdef SINGLE_LINK_DEBUG
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
#endif
    cluster *cnew;
#ifdef SINGLE_LINK_DEBUG
         Rcpp::Rcout << "j1p:" << clust(c1p)
                     << ", j2p:" << clust(c2p)
                << std::endl;
#endif
    if (i == c1->min_d) {
      // we don't need to create a new cluster for c1, we can just modify this
#ifdef SINGLE_LINK_DEBUG
            Rcpp::Rcout << "modifying first cluster: " << clust(c1) << std::endl;
#endif
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
#ifdef SINGLE_LINK_DEBUG
            Rcpp::Rcout << "modifying second cluster: " << clust(c2) << std::endl;
#endif
      cnew = c2;
      remove_child(c1p, c1);
      add_child(cnew, c1);
      c1->parent = cnew;
    } else if (c1p && c1p == c2p && c1p->n_child == 2) {
#ifdef SINGLE_LINK_DEBUG
            Rcpp::Rcout << "modifying common parent cluster: " << clust(c1p) << std::endl;
#endif
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
#ifdef SINGLE_LINK_DEBUG
            Rcpp::Rcout << "finished initializing cluster " << clust(cnew) << std::endl;
#endif
    }
    if (c1p && (c2p == nullptr || c2p->min_d >= c1p->min_d)) {
#ifdef SINGLE_LINK_DEBUG
            Rcpp::Rcout << "parent of new/modified cluster " << clust(cnew)
                        << " should be " << clust(c1p)
                        << std::endl;
#endif
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
#ifdef SINGLE_LINK_DEBUG
            Rcpp::Rcout << "parent of new/modified cluster " << clust(cnew)
                        << " should be " << clust(c2p)
                        << std::endl;
#endif
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
    }
    if (min_d1 < min_d2) i = min_d2; else i = min_d1;
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
#ifdef SINGLE_LINK_TEST
  validate();
#endif
};

void ClusterTree::merge_children(cluster *cdest, cluster *csrc) {
#ifdef SINGLE_LINK_DEBUG
   Rcpp::Rcout << "-merging children of cluster " << clust(csrc)
               << " into cluster " << clust(cdest)
               << std::endl;
#endif
  // CASE 1: the source is null
  // nothing to do
  if (csrc == nullptr) {
#ifdef SINGLE_LINK_DEBUG
         Rcpp::Rcout << " -finished merging children of cluster " << clust(csrc)
                     << " into cluster " << clust(cdest) << " (no-op)" << std::endl;
#endif
    return;
  }
#ifdef SINGLE_LINK_DEBUG
      Rcpp::Rcout << " - source cluster " << clust(csrc)
                  << ": n_child=" << csrc->n_child
                  << ", first_child=" << clust(csrc->first_child)
                  << ", last_child=" << clust(csrc->last_child)
                  << std::endl;
#endif
  // Reassign the parents
  // this is the slow part O(n)
  auto src_child = csrc->first_child;
  if (cdest == nullptr) {
#ifdef SINGLE_LINK_DEBUG
         Rcpp::Rcout << " - destination cluster " << clust(cdest)
                     << " is null" << std::endl
                     << " - removing children from source cluster..." << std::endl;
#endif
    while (src_child) {
#ifdef SINGLE_LINK_FULL_DEBUG
            Rcpp::Rcout << "  - removing child " << clust(src_child) << std::endl;
#endif
      src_child->parent = nullptr;
      src_child->prev_sib = nullptr;
      src_child->next_sib = nullptr;
      src_child = src_child->next_sib;
      csrc->n_child--;
    }
    csrc->first_child = nullptr;
    csrc->last_child = nullptr;
#ifdef SINGLE_LINK_DEBUG
         Rcpp::Rcout << "  -finished removing children."
                     << std::endl
                     << "  - source cluster " << clust(csrc)
                     << ": n_child=" << csrc->n_child
                     << ", first_child=" << clust(csrc->first_child)
                     << ", last_child=" << clust(csrc->last_child)
                     << std::endl;
#endif
  } else {
#ifdef SINGLE_LINK_DEBUG
         Rcpp::Rcout << "  - destination cluster " << clust(cdest)
                     << ": n_child=" << cdest->n_child
                     << ", first_child=" << clust(cdest->first_child)
                     << ", last_child=" << clust(cdest->last_child)
                     << std::endl << " - merging..." << std::endl;
#endif
    // both clusters exist, so do the splice
    auto dest_child = cdest->last_child;
    if (dest_child && src_child) {
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
#ifdef SINGLE_LINK_FULL_DEBUG
            Rcpp::Rcout << "  - transferring child " << clust(src_child) << std::endl;
#endif
      src_child->parent = cdest;
      csrc->n_child--;
      cdest->n_child++;
      src_child = src_child->next_sib;
    }
    csrc->first_child = nullptr;
    csrc->last_child = nullptr;
#ifdef SINGLE_LINK_DEBUG
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
#endif
  }
  return;
};

void ClusterTree::remove_child(cluster *parent, cluster *child) {
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << "-removing child " << clust(child)
                  << " from cluster " << clust(parent)
                  << std::endl;
#endif
  if (parent == nullptr) {
#ifdef SINGLE_LINK_FULL_DEBUG
         Rcpp::Rcout << " -finished removing child " << clust(child)
                     << " from cluster " << clust(parent)
                     << " (no-op)" << std::endl;
#endif
    return;
  }
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << "  -parent " << clust(parent)
                  << ": n_child=" << parent->n_child
                  << ", first_child=" << clust(parent->first_child)
                  << ", last_child=" << clust(parent->last_child)
                  << std::endl
                  << "  -child " << clust(child)
                  << ": prev_sib=" << clust(child->prev_sib)
                  << ", next_sib=" << clust(child->next_sib)
                  << std::endl;
#endif

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
#ifdef SINGLE_LINK_FULL_DEBUG
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
#endif
  return;
};

void ClusterTree::add_child(cluster * parent, cluster * child) {
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << "-adding child " << clust(child)
                  << " to cluster " << clust(parent)
                  << std::endl;
#endif
  if (parent == nullptr) {
#ifdef SINGLE_LINK_FULL_DEBUG
         Rcpp::Rcout << " -finished adding child " << clust(child)
                     << " to cluster " << clust(parent)
                     << " (no-op)" << std::endl;
#endif
    return;
  }
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << "  -parent " << clust(parent)
                  << ": n_child=" << parent->n_child
                  << ", first_child=" << clust(parent->first_child)
                  << ", last_child=" << clust(parent->last_child)
                  << std::endl
                  << "  -child " << clust(child)
                  << ": prev_sib=" << clust(child->prev_sib)
                  << ", next_sib=" << clust(child->next_sib)
                  << std::endl;
#endif
  if (parent->last_child == nullptr) {
    parent->first_child = parent->last_child = child;
    parent->n_child++;
#ifdef SINGLE_LINK_FULL_DEBUG
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
#endif
    return;
  }

  child->prev_sib = parent->last_child;
  parent->last_child->next_sib = child;
  child->next_sib = nullptr;
  parent->last_child = child;
  parent->n_child++;
#ifdef SINGLE_LINK_FULL_DEBUG
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
#endif
  return;
}

void ClusterTree::shift_to_parent(cluster *& c, cluster *& cp) const {
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << "-shifting from child " << clust(c)
                  << " to parent " << clust(cp)
                  << std::endl;
#endif
  c = cp;
  if (cp == nullptr) {
#ifdef SINGLE_LINK_FULL_DEBUG
         Rcpp::Rcout << " -finished shifting to parent " << clust(c)
                     << " (NO_CLUST)"
                     << std::endl;
#endif
    return;
  }
  cp = cp->parent;
#ifdef SINGLE_LINK_FULL_DEBUG
      Rcpp::Rcout << " -finished shifting to parent " << clust(c) << std::endl;
#endif
  return;
}

#ifdef OPTIMOTU_R
int ClusterTree::hclust_ordering(cluster * top, int start, Rcpp::IntegerVector &order) const {
  // Rcpp::Rcout << "ordering cluster " << top << " with start=" << start << std::endl;
  int i = start;
  cluster *c = top->first_child;
  if (c == nullptr) {
    // Rcpp::Rcout << i << ": " << top << std::endl;
    order[i++] = (top - this->my_pool) + 1;
  } else {
    do {
      // Rcpp::Rcout << "getting cluster " << j << std::endl;
      if (c < my_pool + n) {
        // Rcpp::Rcout << i << ": " << j << std::endl;
        order[i++] = c - my_pool + 1;
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
#endif

void ClusterTree::merge_into(DistanceConsumer &consumer) const {
  // std::lock_guard<std::mutex> lock(this->mutex);
  tbb::queuing_rw_mutex::scoped_lock lock(this->mutex, true);
  for (auto c = this->my_pool + this->n; c < this->my_pool + 2*this->n; ++c) {
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

void ClusterTree::merge_into(ClusterAlgorithm &consumer) const {
  // std::lock_guard<std::mutex> lock{this->mutex};
  tbb::queuing_rw_mutex::scoped_lock lock(this->mutex, true);
  for (auto c = this->my_pool + this->n; c < this->my_pool + 2*this->n; ++c) {
    if (c->allocated) {
      if (c->first_child) {
        cluster * next = c->first_child->next_sib;
        while (next) {
          consumer(c->first_child->id, next->id, c->min_d);
          next = next->next_sib;
        }
      }
    }
  }
}

double ClusterTree::max_relevant(j_t seq1, j_t seq2) const {
  // std::lock_guard<std::mutex> lock(this->mutex);
  tbb::queuing_rw_mutex::scoped_lock lock(this->mutex, true);
  cluster* c1 = this->get_cluster(seq1);
  cluster* c1p = c1->parent;
  cluster* c2 = this->get_cluster(seq2);
  cluster* c2p = c2->parent;
  // search from the starting tips up the tree, until either the clusters
  // match (in which case nothing will happen) or both clusters are at the
  // correct height.
  d_t max_d1 = c1->max_d();
  d_t max_d2 = c2->max_d();
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



#ifdef OPTIMOTU_R
void ClusterTree::write_to_matrix(RcppParallel::RMatrix<int> &out) {
  // std::lock_guard<std::mutex> lock(this->mutex);
  tbb::queuing_rw_mutex::scoped_lock lock(this->mutex, true);
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

Rcpp::List ClusterTree::as_hclust(const Rcpp::CharacterVector &seqnames) const {
  // std::lock_guard<std::mutex> lock{this->mutex};
  tbb::queuing_rw_mutex::scoped_lock lock(this->mutex, true);
  Rcpp::IntegerMatrix merge(this->n - 1, 2);
  Rcpp::NumericVector height(this->n-1);
  Rcpp::IntegerVector order(this->n);

  // Rcpp::Rcout << "making cluster/distance pairs" << std::endl;
  // depths and indices of the valid clusters
  std::vector<std::pair<d_t, cluster*>> cluster_dj;
  for (cluster *c = this->my_pool; c < this->my_pool + 2*this->n; ++c) {
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
      if (left - this->my_pool < this->n && first_pass) {
        // Rcpp::Rcout << "left child = " << -left << std::endl;
        *merge1 = this->my_pool - left - 1;
      } else if (first_pass) {
        // Rcpp::Rcout << "left child = " << row_key[child_j] << std::endl;
        *merge1 = row_key[child_c];
      } else {
        // Rcpp::Rcout << "left child = " << row << std::endl;
        *merge1 = row; //not incremented yet; this is the previous row
      }
      if (next_c - this->my_pool < this->n) {
        // Rcpp::Rcout << "right child = " << -(int)next_j << std::endl;
        *merge2 = (int)(this->my_pool - next_c - 1);
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
  for (cluster * c = this->my_pool; c < this->my_pool + this->n; c++) {
    if (c->parent == nullptr) {
      parentless.push_back(c);
    }
  }
  if (parentless.size() > 1) {
    // Rcpp::Rcout << "adding dummy parent for parentless clusters" << std::endl;
    int left;
    if (parentless[0] - this->my_pool < this->n) {
      left = (int)(this->my_pool - parentless[0] - 1);
    } else {
      left = row_key[parentless[0]];
    }
    for (size_t i = 1; i < parentless.size(); i++) {
      int right;
      if (parentless[i] - this->my_pool < this->n) {
        right = (int)(this->my_pool - parentless[i] - 1);
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
#endif

ClusterTree * ClusterTree::make_child() {
  // std::lock_guard<std::mutex> lock(this->mutex);
  tbb::queuing_rw_mutex::scoped_lock lock(this->mutex);
  if (own_child) {
    ClusterTree * child = new ClusterTree(this);
    this->children.insert(child);
    return child;
  }
  this->own_child = true;
  return this;
};

#ifdef SINGLE_LINK_TEST
void ClusterTree::validate() const {
  bool err = false;
#ifdef SINGLE_LINK_DEBUG
  Rcpp::Rcerr << "validating..." << std::endl;
#endif
  for (cluster * c = my_pool; c < my_pool + 2*n; ++c) {
    if (!c->allocated) continue;

    if (c->min_d < -1 || (c->min_d >= m && c->min_d != NO_DIST)) {
      Rcpp::Rcerr << "validation error: cluster " << clust(c)
                  << " has invalid min_d: " << c->min_d << std::endl;
      err = true;
    }
    std::int32_t max = c->max_d();
    if (max < 0 || (max >= m && max != NO_DIST)) {
      Rcpp::Rcerr << "validation error: cluster " << clust(c)
                  << " has invalid max_d: " << max << std::endl;
      err = true;
    }
    if (max < c->min_d) {
      Rcpp::Rcerr << "validation error: cluster " << clust(c)
                  << " has invalid max_d: " << max
                  << " which is less than its min_d: " << c->min_d << std::endl;
      err = true;
    }
    if (c->parent && (c->parent < my_pool + n || c->parent > my_pool + 2*n)) {
      Rcpp::Rcerr << "validation error: cluster " << clust(c)
                  << " has invalid parent: " << clust(c->parent)
                  << std::endl;
      err = true;
    } else if (c->parent && !c->parent->allocated) {
      Rcpp::Rcerr << "validation error: cluster " << clust(c)
                  << "'s parent " << clust(c->parent) << " is not allocated "
                  << std::endl;
      err = true;
    }
    if (c->n_child == 1) {
      Rcpp::Rcerr << "validation error: cluster " << clust(c)
                  << " has exactly one child." << std::endl;
      err = true;
    }
    if (c -> first_child && (c->first_child < my_pool || c->first_child > my_pool + 2*n)) {
      Rcpp::Rcerr << "validation error: cluster " << clust(c)
                  << " has invalid first_child: " << clust(c->first_child) << std::endl;
      err = true;
    }
    auto next = c->first_child;
    std::uint32_t kids = 0;
    while (next) {
      kids++;
      if (next->parent != c) {
        Rcpp::Rcerr << "validation error: cluster " << clust(c)
                    << "'s child " << clust(next)
                    << " instead has parent " << clust(next->parent)
                    << std::endl;
        err = true;
      }
      if (!next->allocated) {
        Rcpp::Rcerr << "validation error: cluster " << clust(c)
                    << "'s child " << clust(next)
                    << " is not allocated " << std::endl;
        err = true;
      }
      next = next->next_sib;
    }
    if (kids != c->n_child) {
      Rcpp::Rcerr << "validation error: cluster " << clust(c)
                  << " has " << kids
                  << " children but n_child=" << c->n_child
                  << std::endl;
      err = true;
    }
    if (c -> next_sib && (c->next_sib < my_pool || c->next_sib > my_pool + 2*n)) {
      Rcpp::Rcerr << "validation error: cluster " << clust(c)
                  << " has invalid next_sib: " << clust(c->next_sib)
                  << std::endl;
      err = true;
    }
    if (c->next_sib && c->next_sib->parent != c->parent) {
      Rcpp::Rcerr << "validation error: cluster " << clust(c)
                  << " with parent " << clust(c->parent)
                  << " has sibling " <<  clust(c->next_sib)
                  << " which instead has parent " << clust(c->next_sib->parent)
                  << std::endl;
      err = true;
    }
    if (c->next_sib && !c->next_sib->allocated) {
      Rcpp::Rcerr << "validation error: cluster " << clust(c)
                  << " with parent " << clust(c->parent)
                  << " has sibling " << clust(c->next_sib)
                  << " which is not allocated " << std::endl;
      err = true;
    }
    if (c->next_sib && c->next_sib->prev_sib != c) {
      Rcpp::Rcerr << "validation error: cluster " << c - my_pool
                  << " with parent " << clust(c->parent)
                  << " has sibling " << clust(c->next_sib)
                  << " whose previous sibling is " << clust(c->next_sib->prev_sib)
                  << std::endl;
      err = true;
    }
    if (c->next_sib == c ) {
      Rcpp::Rcerr << "validation error: cluster " << clust(c)
                  << " is its own next sibling" << std::endl;
      err = true;
    }
  }
  if (err) Rcpp::stop("found validation errors");
}
#endif

#ifdef OPTIMOTU_R

Rcpp::RObject single_linkage_tree(
    const std::string file,
    const Rcpp::CharacterVector &seqnames,
    const DistanceConverter &dconv,
    const std::string &output_type,
    const d_t m
) {
  const j_t n = seqnames.size();
  d_t i;

  // keep track of which clusters are free
  // clusters 0 to n-1 are always used (they are the tips)
  // to start with n to 2n-2 are free.
  ClusterTree ct(dconv, n, m);

  std::ifstream infile(file);
  j_t seq1, seq2;
  double dist;
  while(infile >> seq1 >> seq2 >> dist) {
    if (seq1 == seq2) continue;
#ifdef SINGLE_LINK_DEBUG
    Rcpp::Rcout << "seq1: " << seq1 << ", seq2:" << seq2 << ", dist:" << dist << std::endl;
#endif
    i = dconv.convert(dist);
    if (i >= m) continue;
    ct(seq1, seq2, i);
  }

  if (output_type == "matrix") {
    Rcpp::IntegerMatrix out(m, n);
    RcppParallel::RMatrix<int> rm_out(out);
    ct.write_to_matrix(rm_out);
    return out;
  } else if (output_type == "hclust") {
    Rcpp::List out = ct.as_hclust(seqnames);
    return out;
  } else {
    Rcpp::stop("'output_type' must be 'matrix' or 'hclust'");
  }
}

//' @export
 // [[Rcpp::export]]
 Rcpp::RObject single_linkage_tree_uniform(
     const std::string file,
     const Rcpp::CharacterVector &seqnames,
     const float dmin,
     const float dmax,
     const float dstep,
     const std::string output_type
 ) {
   const UniformDistanceConverter dconv(dmin, dstep);
   const int m = (int) ceilf((dmax - dmin)/dstep) + 1;
   return single_linkage_tree(file, seqnames, dconv, output_type, m);
 }

//' @export
 // [[Rcpp::export]]
 Rcpp::RObject single_linkage_tree_array(
     const std::string file,
     const Rcpp::CharacterVector &seqnames,
     const std::vector<double> &thresholds,
     const std::string output_type
 ) {
   const ArrayDistanceConverter dconv(thresholds);
   const int m = thresholds.size();
   return single_linkage_tree(file, seqnames, dconv, output_type, m);
 }

//' @export
 // [[Rcpp::export]]
 Rcpp::RObject single_linkage_tree_cached(
     const std::string file,
     const Rcpp::CharacterVector &seqnames,
     const std::vector<double> &thresholds,
     const double precision,
     const std::string output_type
 ) {
   const CachedDistanceConverter dconv(thresholds, precision);
   const int m = thresholds.size();
   return single_linkage_tree(file, seqnames, dconv, output_type, m);
 }

Rcpp::List single_linkage_mtree(
    const std::string file,
    const Rcpp::CharacterVector &seqnames,
    const DistanceConverter &dconv,
    const std::string &output_type,
    const int m,
    const Rcpp::ListOf<Rcpp::CharacterVector> &preclust,
    const size_t threads=1
) {
  RcppThread::ThreadPool::globalInstance().setNumThreads(threads);
  const auto namestrings = Rcpp::as<std::vector<std::string>>(seqnames);
  const auto pcstrings = Rcpp::as<std::vector<std::vector<std::string>>>(preclust);
  MultipleClusterAlgorithm<ClusterTree> algo(dconv, namestrings, pcstrings, m);

  std::ifstream infile(file);
  j_t seq1, seq2;
  double dist;
  while(infile >> seq1 >> seq2 >> dist) {
    if (seq1 == seq2) continue;
    algo(seq1, seq2, dist);
  }
  if (output_type == "matrix")
    {
    Rcpp::List out = Rcpp::no_init(preclust.size());
    std::vector<RcppParallel::RMatrix<int>> out_mirror;
    for(j_t pc = 0; pc < preclust.size(); pc++) {
      size_t pc_len = preclust[pc].size();
      Rcpp::IntegerMatrix outm = Rcpp::no_init(m, pc_len);
      Rcpp::colnames(outm) = preclust[pc];
      out[pc] = outm;
      out_mirror.emplace_back(outm);
    }
    algo.write_to_matrix(out_mirror);
    return out;
  } else if (output_type == "hclust") {
    return algo.as_hclust();
  } else {
    Rcpp::stop("'output_type' must be 'matrix' or 'hclust'");
  }
}

//' @export
 // [[Rcpp::export]]
 Rcpp::List single_linkage_mtree_uniform(
     const std::string file,
     const Rcpp::CharacterVector &seqnames,
     const std::string &output_type,
     const double dmin,
     const double dmax,
     const double dstep,
     const Rcpp::ListOf<Rcpp::CharacterVector> &preclust,
     const size_t threads=1
 ) {
   const UniformDistanceConverter dconv(dmin, dstep);
   const int m = (int) ceilf((dmax - dmin)/dstep) + 1;
   return single_linkage_mtree(file, seqnames, dconv, output_type, m, preclust, threads);
 }

//' @export
 // [[Rcpp::export]]
 Rcpp::List single_linkage_mtree_array(
     const std::string file,
     const Rcpp::CharacterVector &seqnames,
     const std::string &output_type,
     const std::vector<double> thresholds,
     const Rcpp::ListOf<Rcpp::CharacterVector> &preclust,
     const size_t threads=1
 ) {
   const ArrayDistanceConverter dconv(thresholds);
   const int m = thresholds.size();
   return single_linkage_mtree(file, seqnames, dconv, output_type, m, preclust, threads);
 }

//' @export
 // [[Rcpp::export]]
 Rcpp::List single_linkage_mtree_cached(
     const std::string file,
     const Rcpp::CharacterVector &seqnames,
     const std::string &output_type,
     const std::vector<double> thresholds,
     const double precision,
     const Rcpp::ListOf<Rcpp::CharacterVector> &preclust,
     const size_t threads=1
 ) {
   const CachedDistanceConverter dconv(thresholds, precision);
   const int m = thresholds.size();
   return single_linkage_mtree(file, seqnames, dconv, output_type, m, preclust, threads);
 }

#endif
