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

// [[Rcpp::export]]
int intersect_length(const std::vector<int> &c, const std::vector<int> &k) {
  return intersect_length_(c.begin(), c.end(), k.begin(), k.end());
}

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
Rcpp::NumericVector fmeasure(
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
