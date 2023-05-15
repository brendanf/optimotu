#include <Rcpp.h>
#include <RcppParallel.h>

struct SizeAndEntropy {
  double H;
  size_t n;
};

inline void initialize_c_counts(
    const Rcpp::IntegerVector c,
    std::vector<std::pair<int, size_t>> &c_sort,
    std::unordered_map<int, SizeAndEntropy> &c_count,
    size_t N
) {
  double Ninv = 1.0/N;
  c_sort.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    c_sort.emplace_back(c[i], i);
    auto counti = c_count.find(c[i]);
    if (counti == c_count.end()) {
      c_count.emplace(c[i], SizeAndEntropy{Ninv, 1});
    } else {
      ++counti->second.n;
      counti->second.H = counti->second.n*Ninv;
    }
  }
  std::sort(c_sort.begin(), c_sort.end());
}

struct MutualInformationWorker : public RcppParallel::Worker
{
  const RcppParallel::RMatrix<int> k;
  const std::vector<std::pair<int, size_t>> &c_sort;
  const std::unordered_map<int, SizeAndEntropy> &c_count;
  const double N;
  RcppParallel::RVector<double> result;

  MutualInformationWorker(
    Rcpp::IntegerMatrix k,
    std::vector<std::pair<int, size_t>> &c_sort,
    std::unordered_map<int, SizeAndEntropy> &c_count,
    Rcpp::NumericVector result
  )
    : k(k), c_sort(c_sort), c_count(c_count), N(k.ncol()), result(result) {};

  void operator()(std::size_t begin, std::size_t end) {
    int c_clust, k_clust;
    double Hij;
    std::unordered_map<int, SizeAndEntropy> k_count;
    std::unordered_map<int, size_t> intersects;
    for (std::size_t j = begin; j < end; ++j) {
      auto kj = k.row(j);
      k_count.clear();
      for (size_t i = 0; i < N; ++i) {
        auto counti = k_count.find(kj[i]);
        if (counti == k_count.end()) {
          k_count.emplace(kj[i], SizeAndEntropy{1.0/N, 1});
        } else {
          ++counti->second.n;
        }
      }
      for (auto &ki : k_count) {
        if (ki.second.n > 1) ki.second.H = ki.second.n/N;
      }
      auto ci = c_sort.begin();
      auto c_end = c_sort.end();
      while(ci < c_end) {
        c_clust = ci->first;
        intersects.clear();
        while(ci < c_end && ci->first == c_clust) {
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
          Hij = (double)counti.second / N;
          result[j] += Hij * log(Hij / c_count.at(c_clust).H / k_count[kj[counti.first]].H);
        }
      }
    }
  }
};

//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mutual_information(
    const Rcpp::IntegerMatrix k,
    const Rcpp::IntegerVector c,
    int threads = 1L
) {
  size_t N = c.size(), m = k.nrow();
  if (N != k.ncol())
    Rcpp::stop("test clusters 'k' (%d) and true clusters 'c' (%d) must have"
                 " the same number of objects.", k.ncol(), N);
  Rcpp::NumericVector mi(m, 0.0);
  std::vector<std::pair<int, size_t>> c_sort;
  std::unordered_map<int, SizeAndEntropy> c_count;
  initialize_c_counts(c, c_sort, c_count, N);

  MutualInformationWorker worker(k, c_sort, c_count, mi);
  if (threads == 1L) {
    worker(0, m);
  } else {
    RcppParallel::parallelFor(0, m, worker, 1, threads);
  }
  return mi;
}

struct AdjustedMutualInformationWorker : public RcppParallel::Worker
{
  const RcppParallel::RMatrix<int> k;
  const std::vector<std::pair<int, size_t>> &c_sort;
  const std::unordered_map<int, SizeAndEntropy> &c_count;
  const size_t N;
  RcppParallel::RVector<double> mi, ami;
  const std::vector<double> lfact;

  std::vector<double> init_lfact() {
    std::vector<double> v;
    v.reserve(N + 1);
    double s = 0.0;
    v.push_back(s);
    for (size_t i = 1; i <= N; ++i) {
      s += log(i);
      v.push_back(s);
    }
    return v;
  }

  AdjustedMutualInformationWorker(
    Rcpp::IntegerMatrix k,
    std::vector<std::pair<int, size_t>> &c_sort,
    std::unordered_map<int, SizeAndEntropy> &c_count,
    Rcpp::NumericVector mi,
    Rcpp::NumericVector ami
  )
    : k(k), c_sort(c_sort), c_count(c_count), N(k.ncol()), lfact(init_lfact()),
      mi(mi), ami(ami) {};

  void operator()(std::size_t begin, std::size_t end) {
    int c_clust, k_clust;
    double Hij, Hc = 0.0, Hk, Ej, maxH;
    std::unordered_map<int, SizeAndEntropy> k_count;
    std::unordered_map<int, size_t> intersects;
    for (auto ci : c_count) {
      Hc -= ci.second.H * log(ci.second.H);
    }
    for (std::size_t j = begin; j < end; j++) {
      auto kj = k.row(j);
      k_count.clear();
      Hk = 0.0;
      for (size_t i = 0; i < N; ++i) {
        auto counti = k_count.find(kj[i]);
        if (counti == k_count.end()) {
          k_count.emplace(kj[i], SizeAndEntropy{1.0/N, 1});
        } else {
          ++counti->second.n;
        }
      }
      for (auto &ki : k_count) {
        if (ki.second.n > 1) ki.second.H = ki.second.n/N;
        Hk -= ki.second.H * log(ki.second.H);
      }
      maxH = Hk > Hc ? Hk : Hc;
      auto ci = c_sort.begin();
      auto c_end = c_sort.end();
      Ej = 0;
      while(ci < c_end) {
        c_clust = ci->first;
        intersects.clear();
        while(ci < c_end && ci->first == c_clust) {
          k_clust = kj[ci->second];
          auto counti = intersects.find(k_clust);
          if (counti == intersects.end()) {
            intersects.emplace(k_clust, 1);
          } else {
            ++counti->second;
          }
          ++ci;
        }
        for (const auto &counti : intersects) {
          Hij = (double)counti.second / N;
          mi[j] += Hij * log(Hij / c_count.at(c_clust).H / k_count[kj[counti.first]].H);
        }
        size_t ai = c_count.at(c_clust).n;
        double apart = lgamma(ai + 1.0) + lgamma(N - ai + 1.0) - lgamma(N + 1.0);
        for (const auto &ki : k_count) {
          size_t bj = ki.second.n;
          size_t nijmin = 1;
          if (ai + bj > N) nijmin = ai + bj - N;
          size_t nijmax = ai > bj ? bj : ai;
          double bpart = apart + lgamma(bj + 1.0) + lgamma(N - bj + 1.0);
          for (size_t nij = nijmin; nij <= nijmax; nij++) {
            Ej += nij/N * log(N * nij / ai / bj) * exp(
              bpart - lgamma(nij + 1.0) - lgamma(ai - nij + 1.0)
            - lgamma(bj - nij + 1.0) - lgamma(N - ai - bj + nij + 1.0)
            );
          }
        }
      }
      ami[j] = (mi[j] - Ej) / (maxH - Ej);
    }
  }
};

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame adjusted_mutual_information(
    const Rcpp::IntegerMatrix k,
    const Rcpp::IntegerVector c,
    int threads = 1L
) {
  size_t N = c.size(), m = k.nrow();
  if (N != k.ncol())
    Rcpp::stop("test clusters 'k' (%d) and true clusters 'c' (%d) must have"
                 " the same number of objects.", k.ncol(), N);
  Rcpp::NumericVector mi(m, 0.0), ami(m, 0.0);
  std::vector<std::pair<int, size_t>> c_sort;
  std::unordered_map<int, SizeAndEntropy> c_count;
  initialize_c_counts(c, c_sort, c_count, N);

  AdjustedMutualInformationWorker worker(k, c_sort, c_count, mi, ami);
  if (threads == 1) {
    worker(0, m);
  } else {
    RcppParallel::parallelFor(0, m, worker, 1, threads);
  }
  auto out = Rcpp::DataFrame::create(
    Rcpp::Named("MI") = mi,
    Rcpp::Named("AMI") = ami
  );
  return out;
}
