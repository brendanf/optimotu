#ifdef OPTIMOTU_R
#include "optimotu.h"
#include <Rcpp.h>
#include <RcppParallel.h>

inline void initialize_counts(
    const Rcpp::IntegerVector c,
    std::vector<std::pair<int, size_t>> &c_sort,
    std::unordered_map<int, size_t> &c_count,
    size_t N
) {
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
  ) : k(k), c_sort(c_sort), c_count(c_count), tp(tp), fp(fp), fn(fn), tn(tn),
  N(c_sort.size()), Npairs(N * (N-1) / 2) {};

  void operator()(size_t begin, size_t end) {
    int c_clust, k_clust;
    size_t tpj, fpj, fnj;
    std::unordered_map<int, size_t> k_count;
    std::unordered_map<int, size_t> intersects;
    for (std::size_t j = begin; j < end; ++j) {
      auto kj = k.row(j);
      tpj = 0;
      fpj = 0;
      fnj = 0;
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
      while(ci != c_end) {
        c_clust = ci->first;
        intersects.clear();
        while(ci != c_end && ci->first == c_clust) {
          k_clust = kj[ci->second];
          auto counti = intersects.find(k_clust);
          if (counti == intersects.end()) {
            intersects.emplace(k_clust, 1);
          } else {
            ++(counti->second);
          }
          ++ci;
        }
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

//' Confusion matrix for a set of "test" partitions vs. a "true" partition
//'
//' One way to analyze the comparison of two different partitions of the same
//' data is to treat it as a binary classification problem operating on pairs
//' of objects, where pairs should be classified as belonging to the same
//' cluster or not.  Then a "true positive" is a pair which belong to the same
//' cluster in both the "test" partition and the "true" partition; a "false
//' positive" is a pair which belongs to the same cluster in the "test"
//' partition but not the "true" partition; a "false negative" is a pair which
//' belongs to the same cluster in the "true" partition but not in the "test"
//' partition; and a "true" negative is a pair which does not belong to the same
//' cluster in either the "test" partition or the "true" partition.
//'
//' This formulation allows various measures of binary classification
//' performance to be applied to the case of clustering.
//'
//' @return (`data.frame` with `m` rows) for each "test" partition, the number
//' of true positives ("`TP`"), false positives ("`FP`"), false negatives
//' ("`FN`"), and true negatives ("`TN`") relative to the "true" partition.
//'
//' @export
//' @inheritParams mutual_information
// [[Rcpp::export]]
Rcpp::DataFrame confusion_matrix(
    const Rcpp::IntegerMatrix k,
    const Rcpp::IntegerVector c,
    const int threads = 1
) {
  size_t N = c.size(), m = k.nrow();
  if (N != (size_t)k.ncol())
    OPTIMOTU_STOP("test clusters 'k' (%d) and true clusters 'c' (%d) must have"
                 " the same number of objects.", k.ncol(), N);
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
  if (threads == 1) {
    worker(0, m);
  } else {
    RcppParallel::parallelFor(0, m, worker, 1, threads);
  }

  auto out = Rcpp::DataFrame::create(
    Rcpp::Named("TP") = true_positive,
    Rcpp::Named("FP") = false_positive,
    Rcpp::Named("FN") = false_negative,
    Rcpp::Named("TN") = true_negative
  );
  if (k.hasAttribute("dimnames")) {
    Rcpp::List dimnames = k.attr("dimnames");
    Rcpp::RObject rownames = dimnames[0];
    if (rownames != R_NilValue) {
      out.attr("row.names") = Rcpp::as<Rcpp::List>(k.attr("dimnames"))[0];
    }
  }

  return out;
}

#endif // OPTIMOTU_R
