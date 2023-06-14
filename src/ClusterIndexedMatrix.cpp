#include <deque>
#include "ClusterIndexedMatrix.h"

template <class A>
void ClusterIndexedMatrix<A>::initialize() {
  auto c = ca;
  index = new tip[n];
  buffer = new int[m];
  // OPTIMOTU_COUT << "cluster array:" << ca << ":" << ca + m*n - 1
  //           << std::endl
  //           << "index:" << index << ":" << index + n - 1
  //           << std::endl
  //           << "buffer:" << buffer << ":" << buffer + m - 1
  //           << std::endl;
  tip * t = index;
  for (j_t j = 0; j < n; j++) {
    t->j = j;
    t->prev = t - 1; // invalid for j=0; fix below
    t->next = t + 1; // invalid for j=n-1; fix below
    t->column = c;
    ++t;
    for (d_t i = 0; i < m; i++) {
      *c = j;
      ++c;
    }
  }
  index[0].prev = nullptr; // fix from above
  index[n-1].next = nullptr; // fix from above
}

template <>
ClusterIndexedMatrix<>::ClusterIndexedMatrix(SingleClusterAlgorithm * parent) :
  SingleClusterAlgorithm(parent), clust_array(m*n), ca(&clust_array[0]) {
  initialize();
}

template<class A>
ClusterIndexedMatrix<A>::ClusterIndexedMatrix(
    const DistanceConverter &dconv,
    size_t n
) : SingleClusterAlgorithm(dconv, n),
clust_array(dconv.m*n),
ca(&clust_array[0])
{
  initialize();
}

template<class A>
ClusterIndexedMatrix<A>::ClusterIndexedMatrix(
  const DistanceConverter &dconv,
  init_matrix_t &im
) : SingleClusterAlgorithm(dconv, im), clust_array(im), ca(&clust_array[0]) {
  initialize();
}

template <class A>
ClusterIndexedMatrix<A>::~ClusterIndexedMatrix() {
  delete[] index;
  delete[] buffer;
}

template <class A>
void ClusterIndexedMatrix<A>::dump_index() {
  OPTIMOTU_CERR << std::endl << "index dump:" << std::endl;
  for (int i = 0; i < n; i++) {
    OPTIMOTU_CERR << "tip=" << i << " prev=";
    if (index[i].prev) {
      OPTIMOTU_CERR << index[i].prev->j << "(";
      if (index[i].prev_d == NO_DIST) {
        OPTIMOTU_CERR << "Inf";
      } else {
        OPTIMOTU_CERR << index[i].prev_d << ")";
      }
    } else {
      OPTIMOTU_CERR << "NA(NA)";
    }
    OPTIMOTU_CERR << " next=";
    if (index[i].next) {
      OPTIMOTU_CERR << index[i].next->j << "(";
      if (index[i].next_d == NO_DIST) {
        OPTIMOTU_CERR << "Inf";
      } else {
        OPTIMOTU_CERR << index[i].next_d << ")";
      }
    } else {
      OPTIMOTU_CERR << "NA(NA)";
    }
    OPTIMOTU_CERR << std::endl;
  }
}

template <class A>
void ClusterIndexedMatrix<A>::print_index() {
  auto t = index;
  OPTIMOTU_COUT << t->j;
  do {
    t=t->next;
    if (t->prev_d == NO_DIST) {
      OPTIMOTU_COUT << "|" << t->j;
    } else {
      OPTIMOTU_COUT << "(" << t->prev_d << ")" << t->j;
    }
  } while (t->next);
  OPTIMOTU_COUT << std::endl;
}

template <class A>
void ClusterIndexedMatrix<A>::verify_index() {
  auto t = index;
  int i = 0;
  while (t->next) {
    t = t->next;
    i++;
    if (i == n) {
      dump_index();
      Rcpp::stop("circular reference formed");
    }
  }
  if (i < n - 1) {
    dump_index();
    Rcpp::stop("broken list");
  }
}

template <class A>
void ClusterIndexedMatrix<A>::heal_splice() {
  d_t join_d = std::max(trev.prev_d, tfwd.next_d);
  if (trev.prev) {
    trev.prev->next = tfwd.next;
    trev.prev->next_d = join_d;
  }
  if (tfwd.next) {
    tfwd.next->prev= trev.prev;
    tfwd.next->prev_d = join_d;
  }
}

template <class A>
bool ClusterIndexedMatrix<A>::index_splice(tip *&t1max, tip *&t2min, tip *&t2max, d_t i) {
  // OPTIMOTU_COUT << "splicing t1max=" << t1max
  //           << " t2min=" << t2min
  //           << " t2max=" << t2max
  //           << std::endl;
  if (t2min == &trev) {
    // OPTIMOTU_COUT << "t2min == trev" << std::endl;
    if (t2max == &tfwd) {
      // OPTIMOTU_COUT << "t2max == tfwd" << std::endl;
      return false;
    } else {
      // OPTIMOTU_COUT << "t2max != tfwd" << std::endl;
      t2min = t2min->prev->next;
    }
  } else if (t2max == &tfwd) {
    // OPTIMOTU_COUT << "t2min != trev" << std::endl;
    // OPTIMOTU_COUT << "t2max == tfwd" << std::endl;
    if (t2max->next) {
      // OPTIMOTU_COUT << "t2max->next: " << t2max->next << std::endl;
      // OPTIMOTU_COUT << "t2max->next->prev: " << t2max->next->prev << std::endl;
      t2max = t2max->next->prev;
    } else {
      t2max = index + n - 1;
    }
  } else if (t1max == &tfwd) {
    // OPTIMOTU_COUT << "t2min != trev" << std::endl;
    // OPTIMOTU_COUT << "t1max == tfwd" << std::endl;
    if (t1max ->next) {
      // OPTIMOTU_COUT << "t1max->next: " << t1max->next << std::endl;
      // OPTIMOTU_COUT << "t1max->next->prev: " << t1max->next->prev << std::endl;
      t1max = t1max->next->prev;
    } else {
      t1max = index + n - 1;
    }
  }
  if (t1max->next == t2min) {
    // OPTIMOTU_COUT << "shortcut merge between " << t1max->j << " and "
    //           << t2min->j << " at " << i << std::endl;
    t1max->next_d = i;
    t2min->prev_d = i;
    t1max = t2max;
    return true;
  }
  // OPTIMOTU_COUT << "Before splice:" << std::endl << "t1max=";
  // if (t1max->j == NO_CLUST) {
  // OPTIMOTU_COUT << "NA";
  // } else {
  // OPTIMOTU_COUT << t1max->j;
  // }
  // OPTIMOTU_COUT << " t1next=";
  // if (t1max->next == nullptr || t1max->next->j == NO_CLUST) {
  // OPTIMOTU_COUT << "NA";
  // } else {
  // OPTIMOTU_COUT << t1max->next->j;
  // }
  // OPTIMOTU_COUT << " t2min=";
  // if (t2min->j == NO_CLUST) {
  // OPTIMOTU_COUT << "NA";
  // } else {
  // OPTIMOTU_COUT << t2min->j;
  // }
  // OPTIMOTU_COUT << " t2prev=";
  // if (t2min->prev == nullptr || t2min->prev->j == NO_CLUST) {
  // OPTIMOTU_COUT << "NA";
  // } else {
  // OPTIMOTU_COUT << t2min->prev->j;
  // }
  // OPTIMOTU_COUT << " t2max=";
  // if (t2max->j == NO_CLUST) {
  // OPTIMOTU_COUT << "NA";
  // } else {
  // OPTIMOTU_COUT << t2max->j;
  // }
  // OPTIMOTU_COUT << " t2next=";
  // if (t2max->next == nullptr || t2max->next->j == NO_CLUST) {
  // OPTIMOTU_COUT << "NA";
  // } else {
  // OPTIMOTU_COUT << t2max->next->j;
  // }
  // OPTIMOTU_COUT << std::endl;
  tfwd.next   = t2max->next;
  tfwd.next_d = t2max->next_d;
  trev.prev   = t2min->prev;
  trev.prev_d = t2min->prev_d;
  t2max->next   = t1max->next;
  t2max->next_d = t1max->next_d;
  if (t1max->next) {
    t1max->next->prev= t2max;
  }
  t1max->next   = t2min;
  t1max->next_d = i;
  t2min->prev   = t1max;
  t2min->prev_d = i;
  t1max = t2max;
  t2max = &tfwd;
  t2min = &trev;
  heal_splice();
  // OPTIMOTU_COUT << "After splice:" << std::endl << "t1max=";
  // if (t1max->j == NO_CLUST) {
  // OPTIMOTU_COUT << "NA";
  // } else {
  // OPTIMOTU_COUT << t1max->j;
  // }
  // OPTIMOTU_COUT << " t1next=";
  // if (t1max->next == nullptr || t1max->next->j == NO_CLUST) {
  // OPTIMOTU_COUT << "NA";
  // } else {
  // OPTIMOTU_COUT << t1max->next->j;
  // }
  // OPTIMOTU_COUT << " t2min=";
  // if (t2min->j == NO_CLUST) {
  // OPTIMOTU_COUT << "NA";
  // } else {
  // OPTIMOTU_COUT << t2min->j;
  // }
  // OPTIMOTU_COUT << " t2prev=";
  // if (t2min->prev == nullptr || t2min->prev->j == NO_CLUST) {
  // OPTIMOTU_COUT << "NA";
  // } else {
  //   OPTIMOTU_COUT << t2min->prev->j;
  // }
  // OPTIMOTU_COUT << " t2max=";
  // if (t2max->j == NO_CLUST) {
  //   OPTIMOTU_COUT << "NA";
  // } else {
  //   OPTIMOTU_COUT << t2max->j;
  // }
  // OPTIMOTU_COUT << " t2next=";
  // if (t2max->next == nullptr || t2max->next->j == NO_CLUST) {
  //   OPTIMOTU_COUT << "NA";
  // } else {
  //   OPTIMOTU_COUT << t2max->next->j;
  // }
  // OPTIMOTU_COUT << std::endl;
  return false;
}

template <class A>
void ClusterIndexedMatrix<A>::operator()(j_t seq1, j_t seq2, d_t i, int thread) {
  if (i >= m) return;
  tbb::queuing_rw_mutex::scoped_lock lock{mutex, true};
  if (clust_array[i + seq1*m] == clust_array[i + seq2*m]) {
    return;
  }
  lock.upgrade_to_writer();
  // OPTIMOTU_COUT << std::endl << "#### seq1=" << seq1
  //           << " seq2=" << seq2
  //           << " i=" << i << std::endl;
  if (clust_array[i + seq1*m] == clust_array[i + seq2*m]) {
    // OPTIMOTU_COUT << "no-op, already clustered" << std::endl;
    return;
  }
  // imin == last row index where seq1 and seq2 are in separate clusters
  // imax == first row index where seq1 and seq2 are in the same cluster
  d_t imin = i, imax = m;
  // shortcut for case when the clusters are not yet joined at any level
  if (clust_array[m-1 + seq1*m] != clust_array[m-1 + seq2*m]) {
    imin = m-1;
  }
  int *c1min = ca + imin + seq1*m;
  int *c2min = ca + imin + seq2*m;
  int *c1max = ca + imax + seq1*m;
  // binary search to find imax and imin
  while (c1max - c1min > 1) {
    int diff = (c1max-c1min) / 2;
    int *c1mean = c1min + diff;
    int *c2mean = c2min + diff;
    if (*c1mean == *c2mean) {
      c1max=c1mean;
    } else {
      c1min=c1mean;
      c2min=c2mean;
    }
  }
  imax=c1max - (ca + seq1*m);
  imin=imax-1;
  // OPTIMOTU_COUT << " imax=" << imax << std::endl;
  size_t j1m, j2m;
  if (seq1 > seq2) {
    j1m = seq2*m;
    j2m = seq1*m;
  } else {
    j1m = seq1*m;
    j2m = seq2*m;
  }
  // OPTIMOTU_COUT << "imax=" << imax << std::endl;
  // cache the new clusters for seq1 and seq2 at each level where they
  // need to be joined
  // OPTIMOTU_COUT << "newclust=";
  for (d_t ii = i; ii < imax; ii++) {
    buffer[ii] = std::min(clust_array[ii + j1m], clust_array[ii + j2m]);
    // OPTIMOTU_COUT << buffer[ii] << " ";
  }
  // OPTIMOTU_COUT << "buffer filled" << std::endl;
  // Rcpp::Rcout << "done" << std::endl;
  // for each sequence which is clustered with seq1 or seq2 at imin,
  // assign it to newclust in the same range where it is
  tip *t1min = index + seq1, *t1max = index + seq1, *t2min = index + seq2,
    *t2max = index + seq2;
  // OPTIMOTU_COUT << "filling t1; t1=" << t1min << "=index[" << t1min - index << "]"
  //           << std::endl;
  // OPTIMOTU_COUT << "t1->column=" << t1min->column << "=clust_array["
  //           << (t1min->column - ca) / m << ", " << (t1min->column-ca) % m << "]"
  //           << std::endl;
  memcpy(t1min->column + i, buffer + i, (imax - i) * sizeof(int));
  // OPTIMOTU_COUT << "filling t2; t2=" << t2min << "=index[" << t2min - index << "]"
  //           << std::endl;
  memcpy(t2min->column + i, buffer + i, (imax - i) * sizeof(int));
  bool merged = false;
  while (i < imax) {
    // OPTIMOTU_COUT << "filling at i=" << i << std::endl;
    d_t nexti = imax;
    size_t total = 0;
    while (t1min->prev_d <= i) {
      t1min = t1min->prev;
      // OPTIMOTU_COUT << "filling t1min; t1min=" << t1min << "=index[" << t1min - index << "]"
      //           << std::endl
      //           << buffer + i << ":" << buffer + imax
      //           << " -> " << t1min->column + i << ":" << t1min->column + imax
      //           << std::endl;
      memcpy(t1min->column + i, buffer + i, (imax - i) * sizeof(int));
      if (total++ > n) {
        // OPTIMOTU_COUT << "infinite loop" << std::endl;
        Rcpp::stop("infinite loop");
      }
    }
    // OPTIMOTU_COUT << "finished t1min" << std::endl;
    total = 0;
    nexti = t1min->prev_d < nexti ? t1min->prev_d : nexti;
    while (t1max->next_d <= i) {
      t1max = t1max->next;
      // OPTIMOTU_COUT << "filling t1max; t1max=" << t1max << "=index[" << t1max - index << "]"
      //           << std::endl
      //           << buffer + i << ":" << buffer + imax
      //           << " -> " << t1max->column + i << ":" << t1max->column + imax
      //           << std::endl;
      memcpy(t1max->column + i, buffer + i, (imax - i) * sizeof(int));
      if (total++ > n) {
        // OPTIMOTU_COUT << "infinite loop" << std::endl;
        Rcpp::stop("infinite loop");
      }
    }
    // OPTIMOTU_COUT << "finished t1max" << std::endl;
    total = 0;
    nexti = t1max->next_d < nexti ? t1max->next_d : nexti;
    if (merged) {
      i = nexti;
      continue;
    }

    while (t2min->prev_d <= i) {
      t2min = t2min->prev;
      // OPTIMOTU_COUT << "filling t2min; t2min=" << t2min << "=index[" << t2min - index << "]"
      //           << std::endl
      //           << buffer + i << ":" << buffer + imax
      //           << " -> " << t2min->column + i << ":" << t2min->column + imax
      //           << std::endl;
      memcpy(t2min->column + i, buffer + i, (imax - i) * sizeof(int));
      if (total++ > n) Rcpp::stop("infinite loop");
    }
    // OPTIMOTU_COUT << "finished t2min" << std::endl;
    total = 0;
    nexti = t2min->prev_d < nexti ? t2min->prev_d : nexti;
    while (t2max->next_d <= i) {
      t2max = t2max->next;
      // OPTIMOTU_COUT << "filling t2max; t2max=" << t2max << "=index[" << t2max - index << "]"
      //           << std::endl
      //           << buffer + i << ":" << buffer + imax
      //           << " -> " << t2max->column + i << ":" << t2max->column + imax
      //           << std::endl;
      memcpy(t2max->column + i, buffer + i, (imax - i) * sizeof(int));
      if (total++ > n) Rcpp::stop("infinite loop");
    }
    // OPTIMOTU_COUT << "finished t2max" << std::endl;
    nexti = t2max->next_d < nexti ? t2max->next_d : nexti;

    if (t1min == t2min) {
      // OPTIMOTU_COUT << "no-op, already merged" << std::endl;
      merged = true;
      continue;
    }
    if (t1min->j < t2min->j) {
      merged = index_splice(t1max, t2min, t2max, i);
    } else {
      merged = index_splice(t2max, t1min, t1max, i);
    }
    // OPTIMOTU_COUT << "finished merging" << std::endl;
    i = nexti;
  }
  // OPTIMOTU_COUT << "finished seq1=" << seq1
  //           << " seq2=" << seq2
  //           << " i=" << i << std::endl;
  // verify_index();
  // print_index();
}

template <class A>
void ClusterIndexedMatrix<A>::merge_into(DistanceConsumer &consumer) {
  tbb::queuing_rw_mutex::scoped_lock lock(this->mutex, true);
  tip * t = index;
  while (t->next != nullptr) {
    if (t->next_d != NO_DIST) consumer(t->j, t->next->j, dconv.inverse(t->next_d));
    t = t->next;
  }
}

template <class A>
void ClusterIndexedMatrix<A>::merge_into(ClusterAlgorithm &consumer) {
  tbb::queuing_rw_mutex::scoped_lock lock(this->mutex, true);
  // TODO check that the distance converters are really compatible
  tip * t = index;
  while (t->next != nullptr) {
    if (t->next_d != NO_DIST) consumer(t->j, t->next->j, t->next_d);
    t = t->next;
  }
}

template <class A>
SingleClusterAlgorithm * ClusterIndexedMatrix<A>::make_child(){
  tbb::queuing_rw_mutex::scoped_lock lock(this->mutex);
  if (own_child) {
    auto child_ptr = new ClusterIndexedMatrix<>(this);
    auto child = std::unique_ptr<ClusterAlgorithm>(
      (ClusterAlgorithm*)child_ptr
    );
    this->children.push_back(std::move(child));
    return (SingleClusterAlgorithm *)child_ptr;
  }
  this->own_child = true;
  return this;
}

template <class A>
double ClusterIndexedMatrix<A>::max_relevant(j_t seq1, j_t seq2, int thread) const {
  if (seq1 == seq2) return 0.0;
  tbb::queuing_rw_mutex::scoped_lock lock(mutex, true);
  j_t j1 = seq1*m, j2 = seq2*m;
  auto c1 = ca + j1, c2 = ca + j2,
    c1max = c1+m;
  if (clust_array[j1 + m - 1] != clust_array[j2 + m - 1]) {
    return dconv.inverse(m-1);
  }
  if (*c1 == *c2) return -1.0;
  // linear search version
  while (c1 < c1max && *c1 != *c2) {
    ++c1;
    ++c2;
  }
  // maybe try binary search?
  return dconv.inverse(c1 - ca - j1 - 1);
}

template <class A>
void ClusterIndexedMatrix<A>::write_to_matrix(internal_matrix_t &out) {
  if (intptr_t(&out[0]) == intptr_t(&clust_array[0])) return;
  std::copy(clust_array.begin(), clust_array.end(), out.begin());
}

#ifdef OPTIMOTU_R

struct OrderElement {
  OrderElement * next;
  int i;
};

// This implementation is copied from ClusterMatrix;
// It would probably be possible to do something faster using the
template <class A>
Rcpp::List ClusterIndexedMatrix<A>::as_hclust(
    const Rcpp::CharacterVector &seqnames
) const {
  // output objects
  Rcpp::IntegerMatrix merge(this->n - 1, 2);
  Rcpp::NumericVector height(this->n-1);
  Rcpp::IntegerVector order(this->n);

  // keep track of which cluster the active tips are in at each stage of clustering
  std::vector<int> clust_id;
  // ordering of later tips in the same cluster
  std::vector<OrderElement> ordering;
  // last element of each cluster (forward_list does not know this)
  std::vector<OrderElement*> last;
  clust_id.reserve(this->n);
  ordering.reserve(this->n);
  last.reserve(this->n);
  // which tips are still "active"
  std::array<std::deque<int>, 2> remaining;
  // keep two versions
  int this_remaining = 0, next_remaining = 1;
  for (int i = 0; i < (int)n;) {
    remaining[this_remaining].push_back(i);
    ++i;
    clust_id.push_back(-i);
    ordering.push_back({NULL, i});
    last.push_back(&ordering.back());
  }
  int last_clust = 0;
  {
    tbb::queuing_rw_mutex::scoped_lock lock{mutex, true};
    for (int j = 0; j < this->m * this->n; j += this->n) {
      double d = this->dconv.inverse(j);
      for (int i : remaining[this_remaining]) {
        int clust = this->ca[j+i];
        if (clust == i) {
          remaining[next_remaining].push_back(i);
        } else {
          merge(last_clust, 0) = clust_id[i];
          merge(last_clust, 1) = clust_id[clust];
          height[last_clust] = d;
          last[clust]->next = &ordering[i];
          last[clust] = last[i];
          clust_id[clust] = ++last_clust;
        }
        remaining[this_remaining].pop_front();
      }
      next_remaining++;
      next_remaining %= 2;
      this_remaining++;
      this_remaining %= 2;
    }
  }
  int j = 0;
  for (int i : remaining[this_remaining]) {
    if (i > 0) {
      merge(last_clust, 0) = i;
      merge(last_clust, 1) = clust_id[0];
      height[last_clust] = 1.0;
      clust_id[0] = ++last_clust;
    }
    OrderElement * e = &ordering[i];
    while (e != NULL) {
      order[j] = e->i;
      ++j;
      e = e->next;
    }
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
