#include "ClusterIndexedMatrix.h"

template <class A>
void ClusterIndexedMatrix<A>::initialize() {
  auto c = ca;
  index = new tip[n];
  buffer = new int[m];
  // std::cout << "cluster array:" << ca << ":" << ca + m*n - 1
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
};

template <class A>
ClusterIndexedMatrix<A>::ClusterIndexedMatrix(ClusterAlgorithm * parent) :
  ClusterAlgorithm(parent), clust_array(m*n), ca(&clust_array[0]) {
  initialize();
}

template <class A>
ClusterIndexedMatrix<A>::~ClusterIndexedMatrix() {
  delete[] index;
  delete[] buffer;
}

template <class A>
void ClusterIndexedMatrix<A>::dump_index() {
  std::cerr << std::endl << "index dump:" << std::endl;
  for (int i = 0; i < n; i++) {
    std::cerr << "tip=" << i << " prev=";
    if (index[i].prev) {
      std::cerr << index[i].prev->j << "(";
      if (index[i].prev_d == NO_DIST) {
        std::cerr << "Inf";
      } else {
        std::cerr << index[i].prev_d << ")";
      }
    } else {
      std::cerr << "NA(NA)";
    }
    std::cerr << " next=";
    if (index[i].next) {
      std::cerr << index[i].next->j << "(";
      if (index[i].next_d == NO_DIST) {
        std::cerr << "Inf";
      } else {
        std::cerr << index[i].next_d << ")";
      }
    } else {
      std::cerr << "NA(NA)";
    }
    std::cerr << std::endl;
  }
}

template <class A>
void ClusterIndexedMatrix<A>::print_index() {
  auto t = index;
  std::cout << t->j;
  do {
    t=t->next;
    if (t->prev_d == NO_DIST) {
      std::cout << "|" << t->j;
    } else {
      std::cout << "(" << t->prev_d << ")" << t->j;
    }
  } while (t->next);
  std::cout << std::endl;
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
  // std::cout << "splicing t1max=" << t1max
  //           << " t2min=" << t2min
  //           << " t2max=" << t2max
  //           << std::endl;
  if (t2min == &trev) {
    // std::cout << "t2min == trev" << std::endl;
    if (t2max == &tfwd) {
      // std::cout << "t2max == tfwd" << std::endl;
      return false;
    } else {
      // std::cout << "t2max != tfwd" << std::endl;
      t2min = t2min->prev->next;
    }
  } else if (t2max == &tfwd) {
    // std::cout << "t2min != trev" << std::endl;
    // std::cout << "t2max == tfwd" << std::endl;
    if (t2max->next) {
      // std::cout << "t2max->next: " << t2max->next << std::endl;
      // std::cout << "t2max->next->prev: " << t2max->next->prev << std::endl;
      t2max = t2max->next->prev;
    } else {
      t2max = index + n - 1;
    }
  } else if (t1max == &tfwd) {
    // std::cout << "t2min != trev" << std::endl;
    // std::cout << "t1max == tfwd" << std::endl;
    if (t1max ->next) {
      // std::cout << "t1max->next: " << t1max->next << std::endl;
      // std::cout << "t1max->next->prev: " << t1max->next->prev << std::endl;
      t1max = t1max->next->prev;
    } else {
      t1max = index + n - 1;
    }
  }
  if (t1max->next == t2min) {
    // std::cout << "shortcut merge between " << t1max->j << " and "
    //           << t2min->j << " at " << i << std::endl;
    t1max->next_d = i;
    t2min->prev_d = i;
    t1max = t2max;
    return true;
  }
  // std::cout << "Before splice:" << std::endl << "t1max=";
  // if (t1max->j == NO_CLUST) {
  // std::cout << "NA";
  // } else {
  // std::cout << t1max->j;
  // }
  // std::cout << " t1next=";
  // if (t1max->next == nullptr || t1max->next->j == NO_CLUST) {
  // std::cout << "NA";
  // } else {
  // std::cout << t1max->next->j;
  // }
  // std::cout << " t2min=";
  // if (t2min->j == NO_CLUST) {
  // std::cout << "NA";
  // } else {
  // std::cout << t2min->j;
  // }
  // std::cout << " t2prev=";
  // if (t2min->prev == nullptr || t2min->prev->j == NO_CLUST) {
  // std::cout << "NA";
  // } else {
  // std::cout << t2min->prev->j;
  // }
  // std::cout << " t2max=";
  // if (t2max->j == NO_CLUST) {
  // std::cout << "NA";
  // } else {
  // std::cout << t2max->j;
  // }
  // std::cout << " t2next=";
  // if (t2max->next == nullptr || t2max->next->j == NO_CLUST) {
  // std::cout << "NA";
  // } else {
  // std::cout << t2max->next->j;
  // }
  // std::cout << std::endl;
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
  // std::cout << "After splice:" << std::endl << "t1max=";
  // if (t1max->j == NO_CLUST) {
  // std::cout << "NA";
  // } else {
  // std::cout << t1max->j;
  // }
  // std::cout << " t1next=";
  // if (t1max->next == nullptr || t1max->next->j == NO_CLUST) {
  // std::cout << "NA";
  // } else {
  // std::cout << t1max->next->j;
  // }
  // std::cout << " t2min=";
  // if (t2min->j == NO_CLUST) {
  // std::cout << "NA";
  // } else {
  // std::cout << t2min->j;
  // }
  // std::cout << " t2prev=";
  // if (t2min->prev == nullptr || t2min->prev->j == NO_CLUST) {
  // std::cout << "NA";
  // } else {
  //   std::cout << t2min->prev->j;
  // }
  // std::cout << " t2max=";
  // if (t2max->j == NO_CLUST) {
  //   std::cout << "NA";
  // } else {
  //   std::cout << t2max->j;
  // }
  // std::cout << " t2next=";
  // if (t2max->next == nullptr || t2max->next->j == NO_CLUST) {
  //   std::cout << "NA";
  // } else {
  //   std::cout << t2max->next->j;
  // }
  // std::cout << std::endl;
  return false;
}

template <class A>
void ClusterIndexedMatrix<A>::operator()(j_t seq1, j_t seq2, d_t i) {
  if (i >= m) return;
  tbb::queuing_rw_mutex::scoped_lock lock{mutex, true};
  if (clust_array[i + seq1*m] == clust_array[i + seq2*m]) {
    return;
  }
  lock.upgrade_to_writer();
  // std::cout << std::endl << "#### seq1=" << seq1
  //           << " seq2=" << seq2
  //           << " i=" << i << std::endl;
  if (clust_array[i + seq1*m] == clust_array[i + seq2*m]) {
    // std::cout << "no-op, already clustered" << std::endl;
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
  // std::cout << " imax=" << imax << std::endl;
  size_t j1m, j2m;
  if (seq1 > seq2) {
    j1m = seq2*m;
    j2m = seq1*m;
  } else {
    j1m = seq1*m;
    j2m = seq2*m;
  }
  // std::cout << "imax=" << imax << std::endl;
  // cache the new clusters for seq1 and seq2 at each level where they
  // need to be joined
  // std::cout << "newclust=";
  for (d_t ii = i; ii < imax; ii++) {
    buffer[ii] = std::min(clust_array[ii + j1m], clust_array[ii + j2m]);
    // std::cout << buffer[ii] << " ";
  }
  // std::cout << "buffer filled" << std::endl;
  // Rcpp::Rcout << "done" << std::endl;
  // for each sequence which is clustered with seq1 or seq2 at imin,
  // assign it to newclust in the same range where it is
  tip *t1min = index + seq1, *t1max = index + seq1, *t2min = index + seq2,
    *t2max = index + seq2;
  // std::cout << "filling t1; t1=" << t1min << "=index[" << t1min - index << "]"
  //           << std::endl;
  // std::cout << "t1->column=" << t1min->column << "=clust_array["
  //           << (t1min->column - ca) / m << ", " << (t1min->column-ca) % m << "]"
  //           << std::endl;
  memcpy(t1min->column + i, buffer + i, (imax - i) * sizeof(int));
  // std::cout << "filling t2; t2=" << t2min << "=index[" << t2min - index << "]"
  //           << std::endl;
  memcpy(t2min->column + i, buffer + i, (imax - i) * sizeof(int));
  bool merged = false;
  while (i < imax) {
    // std::cout << "filling at i=" << i << std::endl;
    d_t nexti = imax;
    size_t total = 0;
    while (t1min->prev_d <= i) {
      t1min = t1min->prev;
      // std::cout << "filling t1min; t1min=" << t1min << "=index[" << t1min - index << "]"
      //           << std::endl
      //           << buffer + i << ":" << buffer + imax
      //           << " -> " << t1min->column + i << ":" << t1min->column + imax
      //           << std::endl;
      memcpy(t1min->column + i, buffer + i, (imax - i) * sizeof(int));
      if (total++ > n) {
        // std::cout << "infinite loop" << std::endl;
        Rcpp::stop("infinite loop");
      }
    }
    // std::cout << "finished t1min" << std::endl;
    total = 0;
    nexti = t1min->prev_d < nexti ? t1min->prev_d : nexti;
    while (t1max->next_d <= i) {
      t1max = t1max->next;
      // std::cout << "filling t1max; t1max=" << t1max << "=index[" << t1max - index << "]"
      //           << std::endl
      //           << buffer + i << ":" << buffer + imax
      //           << " -> " << t1max->column + i << ":" << t1max->column + imax
      //           << std::endl;
      memcpy(t1max->column + i, buffer + i, (imax - i) * sizeof(int));
      if (total++ > n) {
        // std::cout << "infinite loop" << std::endl;
        Rcpp::stop("infinite loop");
      }
    }
    // std::cout << "finished t1max" << std::endl;
    total = 0;
    nexti = t1max->next_d < nexti ? t1max->next_d : nexti;
    if (merged) {
      i = nexti;
      continue;
    }

    while (t2min->prev_d <= i) {
      t2min = t2min->prev;
      // std::cout << "filling t2min; t2min=" << t2min << "=index[" << t2min - index << "]"
      //           << std::endl
      //           << buffer + i << ":" << buffer + imax
      //           << " -> " << t2min->column + i << ":" << t2min->column + imax
      //           << std::endl;
      memcpy(t2min->column + i, buffer + i, (imax - i) * sizeof(int));
      if (total++ > n) Rcpp::stop("infinite loop");
    }
    // std::cout << "finished t2min" << std::endl;
    total = 0;
    nexti = t2min->prev_d < nexti ? t2min->prev_d : nexti;
    while (t2max->next_d <= i) {
      t2max = t2max->next;
      // std::cout << "filling t2max; t2max=" << t2max << "=index[" << t2max - index << "]"
      //           << std::endl
      //           << buffer + i << ":" << buffer + imax
      //           << " -> " << t2max->column + i << ":" << t2max->column + imax
      //           << std::endl;
      memcpy(t2max->column + i, buffer + i, (imax - i) * sizeof(int));
      if (total++ > n) Rcpp::stop("infinite loop");
    }
    // std::cout << "finished t2max" << std::endl;
    nexti = t2max->next_d < nexti ? t2max->next_d : nexti;

    if (t1min == t2min) {
      // std::cout << "no-op, already merged" << std::endl;
      merged = true;
      continue;
    }
    if (t1min->j < t2min->j) {
      merged = index_splice(t1max, t2min, t2max, i);
    } else {
      merged = index_splice(t2max, t1min, t1max, i);
    }
    // std::cout << "finished merging" << std::endl;
    i = nexti;
  }
  // std::cout << "finished seq1=" << seq1
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
ClusterAlgorithm * ClusterIndexedMatrix<A>::make_child(){
  tbb::queuing_rw_mutex::scoped_lock lock(this->mutex);
  if (own_child) {
    ClusterAlgorithm * child = new ClusterIndexedMatrix<std::vector<int>>(this);
    this->children.insert(child);
    return child;
  }
  this->own_child = true;
  return this;
}

template <class A>
double ClusterIndexedMatrix<A>::max_relevant(j_t seq1, j_t seq2) const {
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
};

#ifdef OPTIMOTU_R
template <class A>
void ClusterIndexedMatrix<A>::write_to_matrix(RcppParallel::RMatrix<int> &out) {
  std::copy(clust_array.begin(), clust_array.end(), out.begin());
}

template<>
ClusterIndexedMatrix<RcppParallel::RMatrix<int>>::ClusterIndexedMatrix(
    const DistanceConverter &dconv, Rcpp::IntegerMatrix &im
) :
  ClusterAlgorithm(dconv, im.ncol(), im.nrow()), clust_array(im), ca(&clust_array[0])
{
  initialize();
};
#endif
