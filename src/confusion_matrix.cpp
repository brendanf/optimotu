#include <Rcpp.h>
#include <RcppParallel.h>

std::vector<std::pair<int, size_t>> cluster_map(const std::vector<int> & c) {
  std::vector<std::pair<int, size_t>> key;
  key.reserve(c.size());
  for (size_t i = 0; i < c.size(); i++) {
    key.emplace_back(c[i], i);
  }
  std::sort(key.begin(), key.end());
  return key;
}

std::vector<std::pair<int, size_t>> cluster_map(const Rcpp::IntegerMatrix::ConstRow & c) {
  std::vector<std::pair<int, size_t>> key;
  key.reserve(c.size());
  for (size_t i = 0; i < c.size(); i++) {
    key.emplace_back(c[i], i);
  }
  std::sort(key.begin(), key.end());
  return key;
}

std::vector<int> cluster_key(const std::vector<std::pair<int, size_t>> map) {
  std::vector<int> key(map.size(), -1);
  int j;
  for (int i = 0; i < map.size(); i++) {
    j = map[i].first;
    if (key[j] == -1) key[j] = i;
  }
  return(key);
}


//' @export
//[[Rcpp::export]]
Rcpp::DataFrame confusion_matrix(
    const Rcpp::IntegerMatrix &cluster_matrix,
    const std::vector<int> &known_clusters
) {
  size_t n_thresh = cluster_matrix.nrow();
  std::vector<double> true_positive;
  true_positive.reserve(n_thresh);
  std::vector<double> false_positive;
  false_positive.reserve(n_thresh);
  std::vector<double> false_negative;
  false_negative.reserve(n_thresh);
  std::vector<double> true_negative;
  true_negative.reserve(n_thresh);
  const size_t N = known_clusters.size();
  const size_t N_pairs = N * (N-1) / 2;

  auto known_map = cluster_map(known_clusters);

  for (size_t i = 0; i < n_thresh; i++) {
    auto cluster_row = cluster_matrix.row(i);
    auto map_i = cluster_map(cluster_row);
    auto key_i = cluster_key(map_i);
    size_t tp = 0, fp = 0, fn = 0;
    for (size_t j = 0; j < N; j++) {
      int c1 = known_map[j].first;
      size_t s1 = known_map[j].second;
      // check in "true" cluster for true positives and false negatives
      // because
      size_t k = j + 1;
      while(k < N && known_map[k].first == c1) {
        if (cluster_row[s1] == cluster_row[known_map[k].second]) {
          ++tp;
        } else {
          ++fn;
        }
        ++k;
      }
      // check in "test" cluster for false positives
      int c2 = cluster_row[s1];
      k = key_i[c2];
      while (k < N && map_i[k].first == c2) {
        size_t s2 = map_i[k].second;
        if (s2 == s1) break;
        if (known_clusters[s1] != known_clusters[s2]) ++fp;
        ++k;
      }
    }
    true_positive.push_back(tp);
    false_positive.push_back(fp);
    false_negative.push_back(fn);
    true_negative.push_back(N_pairs - tp - fp - fn);
  }

  auto out = Rcpp::DataFrame::create(
    Rcpp::Named("TP") = true_positive,
    Rcpp::Named("FP") = false_positive,
    Rcpp::Named("FN") = false_negative,
    Rcpp::Named("TN") = true_negative
  );

  return out;
}

inline void initialize_counts(
    const Rcpp::IntegerVector c,
    std::vector<std::pair<int, size_t>> &c_sort,
    std::unordered_map<int, size_t> &c_count,
    size_t N
) {
  double Ninv = 1.0/N;
  c_sort.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    c_sort.emplace_back(c[i], i);
    auto counti = c_count.find(c[i]);
    if (counti == c_count.end()) {
      c_count.emplace(c[i], 1);
    } else {
      ++counti->second;
    }
  }
  std::sort(c_sort.begin(), c_sort.end());
}


struct ConfusionMatrixWorker : RcppParallel::Worker {
  const RcppParallel::RMatrix<int> k;
  const std::vector<std::pair<int, size_t>> &c_sort;
  const std::unordered_map<int, size_t> &c_count;
  RcppParallel::RVector<double> tp, fp, fn, tn;
  const size_t N, Npairs;

  ConfusionMatrixWorker(
    Rcpp::IntegerMatrix k,
    std::vector<std::pair<int, size_t>> &c_sort,
    std::unordered_map<int, size_t> &c_count,
    Rcpp::NumericVector tp,
    Rcpp::NumericVector fp,
    Rcpp::NumericVector fn,
    Rcpp::NumericVector tn
  ) : k(k), c_sort(c_sort), c_count(c_count), N(c_sort.size()), Npairs(N * (N-1) / 2),
  tp(tp), fp(fp), fn(fn), tn(tn) {};

  void operator()(size_t begin, size_t end) {
    int c_clust, k_clust;
    size_t tpj, fpj, fnj, tnj;
    std::unordered_map<int, size_t> k_count;
    std::unordered_map<int, size_t> intersects;
    for (std::size_t j = begin; j < end; ++j) {
      auto kj = k.row(j);
      tpj = 0;
      fpj = 0;
      fnj = 0;
      tnj = 0;
      k_count.clear();
      for (size_t i = 0; i < N; ++i) {
        auto counti = k_count.find(kj[i]);
        if (counti == k_count.end()) {
          k_count.emplace(kj[i], 1);
        } else {
          ++counti->second;
        }
      }
      auto ci = c_sort.begin();
      auto c_end = c_sort.end();
      while(ci < c_end) {
        c_clust = ci->first;
        intersects.clear();
        while(ci->first == c_clust) {
          k_clust = kj[ci->second];
          auto counti = intersects.find(k_clust);
          if (counti == intersects.end()) {
            intersects.emplace(k_clust, 1);
          } else {
            ++(counti->second);
          }
          ++ci;
        }
        double m = 0;
        for (const auto counti : intersects) {
          tpj += counti.second * (counti.second - 1);
          fpj += counti.second * (k_count.at(counti.first) - counti.second);
          fnj += counti.second * (c_count.at(c_clust) - counti.second);
        }
      }
      tp[j] = tpj / 2.0;
      fp[j] = fpj / 2.0;
      fn[j] = fnj / 2.0;
      tn[j] = Npairs - tp[j] - fp[j] - fn[j];
    }
  }
};

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame confusion_matrix2(
    const Rcpp::IntegerMatrix k,
    const Rcpp::IntegerVector c,
    const int ncpu = 1
) {
  size_t N = c.size(), m = k.nrow();
  Rcpp::NumericVector true_positive(m);
  Rcpp::NumericVector false_positive(m);
  Rcpp::NumericVector false_negative(m);
  Rcpp::NumericVector true_negative(m);

  std::vector<std::pair<int, size_t>> c_sort;
  std::unordered_map<int, size_t> c_count;
  initialize_counts(c, c_sort, c_count, N);

  ConfusionMatrixWorker worker(k, c_sort, c_count,
                               true_positive, false_positive,
                               false_negative, true_negative);
  if (ncpu == 1) {
    worker(0, m);
  } else {
    RcppParallel::parallelFor(0, m, worker, 1, ncpu);
  }

  auto out = Rcpp::DataFrame::create(
    Rcpp::Named("TP") = true_positive,
    Rcpp::Named("FP") = false_positive,
    Rcpp::Named("FN") = false_negative,
    Rcpp::Named("TN") = true_negative
  );

  return out;
}
