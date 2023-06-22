#ifdef OPTIMOTU_R

#include "optimotu.h"
#include <Rcpp.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppThread)]]
#include <RcppThread.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>


template <class It1, class It2>
int intersect_length_(It1 ci, It1 cend, It2 ki, It2 kend) {
   int l = 0;
   while (true) {
      if (*ci < *ki) {
         if (++ci == cend) return l;
      } else if (*ki < *ci) {
         if (++ki == kend) return l;
      } else {
         ++l;
         if (++ci == cend || ++ki == kend) return l;
      }
   }
   return l;
}
//' Size of the intersection between two sorted sets
//'
//' This implementation is much faster that `length(intersect(c, k))`. However
//' it assumes (without checking!) that the sets are sorted.
//'
//' @param c (sorted `integer` or ` character` vector) the first set to compare
//' @param k (sorted `integer` or ` character` vector) the second set to compare
//'
//' @return (`integer` count) the number of elements which occur in both sets
//' @export
// [[Rcpp::export]]
int intersect_length(const std::vector<int> &c, const std::vector<int> &k) {
  return intersect_length_(c.begin(), c.end(), k.begin(), k.end());
}

//' @export
//' @rdname intersect_length
// [[Rcpp::export]]
int intersect_length_string(const std::vector<std::string> &c, const std::vector<std::string> &k) {
   return intersect_length_(c.begin(), c.end(), k.begin(), k.end());
}


// [[Rcpp::export]]
double inner_fmeasure(
      const std::vector<int> &cj,
      const Rcpp::ListOf<Rcpp::IntegerVector> &kpartition,
      const std::vector<int> &nk
) {
   double nc = cj.size();
   double m = 0;
   auto nki = nk.begin();
   for (auto kj : kpartition) {
      double mm = intersect_length(Rcpp::as<std::vector<int>>(kj), cj) / (nc + *(nki++));
      if (mm > m) m = mm;
   }
   return nc * m;
}

struct FmeasureWorker : public RcppParallel::Worker
{
   std::vector<std::vector<std::vector<int>>> k;
   std::vector<std::vector<int>> c;
   size_t c_tot = 0;
   RcppParallel::RVector<double> result;

   FmeasureWorker(Rcpp::ListOf<Rcpp::ListOf<Rcpp::IntegerVector>> &k,
                   Rcpp::ListOf<Rcpp::IntegerVector> &c,
                  Rcpp::NumericVector &result)
      : result(result) {
      this->k = Rcpp::as<std::vector<std::vector<std::vector<int>>>>(k);
      this->c = Rcpp::as<std::vector<std::vector<int>>>(c);
      for(auto ci : this->c) {
         c_tot += ci.size();
      }
   };

   void operator()(std::size_t begin, std::size_t end) {
      for (std::size_t i = begin; i < end; i++) {
         result[i] = 0.0;
         for (auto cj : c) {
            double m = 0;
            for (auto kj : k[i]) {
               double mm = intersect_length_(
                  cj.begin(),
                  cj.end(),
                  kj.begin(),
                  kj.end()
               ) /
                  (double)(cj.size() + kj.size());
               if (mm > m) m = mm;
            }
            result[i] += cj.size() * m;
         }
         result[i] *= 2.0 / c_tot;
      }
   }
};

// [[Rcpp::export]]
Rcpp::NumericVector fmeasure_list(
      Rcpp::ListOf<Rcpp::ListOf<Rcpp::IntegerVector>> k,
      Rcpp::ListOf<Rcpp::IntegerVector> c,
      size_t ncpu = 1
) {
   auto k_int = Rcpp::as<std::vector<std::vector<std::vector<int>>>>(k);
   auto c_int = Rcpp::as<std::vector<std::vector<int>>>(c);
   Rcpp::NumericVector result(k.size());
   FmeasureWorker worker(k, c, result);
   if (ncpu > 1) {
      RcppParallel::parallelFor(0, k.size(), worker, ncpu);
   } else {
      worker(0, k.size());
   }
   return result;
}

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

struct FMeasureWorker2 : public RcppParallel::Worker
{
  const RcppParallel::RMatrix<int> k;
  const std::vector<std::pair<int, size_t>> &c_sort;
  const std::unordered_map<int, size_t> &c_count;
  const double N;
  RcppParallel::RVector<double> result;

  FMeasureWorker2(
    Rcpp::IntegerMatrix k,
    std::vector<std::pair<int, size_t>> &c_sort,
    std::unordered_map<int, size_t> &c_count,
    Rcpp::NumericVector result
  )
    : k(k), c_sort(c_sort), c_count(c_count), N(k.ncol()), result(result) {};

  void operator()(std::size_t begin, std::size_t end) {
    int c_clust, k_clust;
    std::unordered_map<int, size_t> k_count;
    std::unordered_map<int, size_t> intersects;
    for (std::size_t j = begin; j < end; ++j) {
      auto kj = k.row(j);
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
          double mm = counti.second /
            (double)(c_count.at(c_clust) + k_count.at(counti.first));
          if (mm > m) m = mm;
        }
        result[j] += c_count.at(c_clust) * m;
      }
      result[j] *= 2.0 / N;
    }
  }
};


// [[Rcpp::export]]
Rcpp::NumericVector fmeasure_matrix(
  Rcpp::IntegerMatrix k,
  Rcpp::IntegerVector c,
  size_t ncpu = 1
) {
  size_t n = c.size(), m = k.nrow();
  if (n != (size_t)k.ncol())
    OPTIMOTU_STOP("test clusters 'k' (%d) and true clusters 'c' (%d) must have"
                 " the same number of objects.", k.ncol(), n);
  Rcpp::NumericVector fm(m, 0.0);
  std::vector<std::pair<int, size_t>> c_sort;
  std::unordered_map<int, size_t> c_count;
  initialize_counts(c, c_sort, c_count, n);
  FMeasureWorker2 worker(k, c_sort, c_count, fm);
  if (ncpu == 1) {
    worker(0, m);
  } else {
    RcppParallel::parallelFor(0, k.nrow(), worker, 1, ncpu);
  }
  return fm;
}

#endif // OPTIMOTU_R
