#ifdef OPTIMOTU_R
#include <Rcpp.h>
#include <RcppParallel.h>

#include<mutex>

#define NEGINF (-std::numeric_limits<double>::infinity())

double log_plus(double x, double y) {
  if (x == NEGINF) return y;
  if (y == NEGINF) return x;
  if (x == y) return x + std::log(2);
  if (x > y) return x + std::log1p(std::exp(y - x));
  return y + std::log1p(std::exp(x - y));
}

void log_add_to(double &x, double y) {
  if (x == NEGINF) {
    x = y;
  } else if (y == NEGINF) {
    return;
  } else if (x == y) {
    x += std::log(2);
  } else if (x > y) {
    x += std::log1p(std::exp(y - x));
  } else {
    x = y + std::log1p(std::exp(x - y));
  }
}

double log_minus(double x, double y) {
  if (x == NEGINF) return y;
  if (y == NEGINF) return x;
  if (x == y) return x + std::log(2);
  if (x > y) return x + std::log1p(-std::exp(y - x));
  return y + std::log1p(-std::exp(x - y));
}

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
      // count the number of items in each test cluster
      k_count.clear();
      for (size_t i = 0; i < N; ++i) {
        auto counti = k_count.find(kj[i]);
        if (counti == k_count.end()) {
          k_count.emplace(kj[i], SizeAndEntropy{1.0/N, 1});
        } else {
          ++counti->second.n;
        }
      }
      // calculate entropy based on final counts
      for (auto &ki : k_count) {
        if (ki.second.n > 1) ki.second.H = ki.second.n/N;
        // Rcpp::Rcerr << "k" << ki.first
        //             << " has " << ki.second.n
        //             << " members and H=" << ki.second.H << std::endl;
      }
      // traverse true clusters
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
          // Rcpp::Rcerr << "intersections between k" << counti.first
          //             << " and c" << c_clust
          //             << ": " << counti.second << std::endl;
          Hij = (double)counti.second / N;
          // Rcpp::Rcerr << "mutual information terms: "
          //             << "Hij = " << Hij << std::endl
          //             << "Hi = " << c_count.at(c_clust).H << std::endl
          //             << "Hj = " << k_count.at(counti.first).H << std::endl
          //             << "MIij = " << Hij * log(Hij / c_count.at(c_clust).H / k_count.at(counti.first).H)
          //             << std::endl;
          result[j] += Hij * log(Hij / c_count.at(c_clust).H / k_count.at(counti.first).H);
        }
      }
    }
  }
};

//' Calculate similarity for a set of alternate "test" partitions vs. a "true" partition
//'
//' @param k (`m` x `n` `integer` matrix) `m` alternative "test" partitions; each
//' row gives the cluster assignment for the `n` objects. Objects with the same
//' cluster ID are clustered together. Cluster IDs do not need to be
//' consecutive, and they do not need to correspond between different
//' partitions.
//' @param c (`integer` vector of length `n`) "True" partition of the `n`
//' objects.
//' @param threads (`integer` count) number of parallel threads to use.
//'
//' @return (`numeric` vector of length `m`) The similarity measure between each
//' of the alternative partitions and the "true" partition.
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mutual_information(
    const Rcpp::IntegerMatrix k,
    const Rcpp::IntegerVector c,
    int threads = 1L
) {
  size_t N = c.size(), m = k.nrow();
  if (N != (size_t)k.ncol())
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

struct ClusterCount {
  size_t size, // size of each cluster
  n, //number of clusters
  j; // which j this is a part of

  bool operator<(ClusterCount x) {
    return size < x.size;
  }
};

struct AdjustedMutualInformationWorker1 : public RcppParallel::Worker
{
  const RcppParallel::RMatrix<int> k; // rows are test partitions
  const std::vector<std::pair<int, size_t>> &c_sort;
  const std::unordered_map<int, SizeAndEntropy> &c_count;
  const size_t N;
  RcppParallel::RVector<double> mi, Hmax;
  std::vector<ClusterCount> &all_k_counts;
  const double Hc; // entropy of c partition
  std::mutex mutex;

  double init_Hc() {
    double Hc = 0.0;
    // calculate entropy of true partition
    for (auto ci : c_count) {
      Hc -= ci.second.H * log(ci.second.H);
    }
    return Hc;
  }

AdjustedMutualInformationWorker1(
    const Rcpp::IntegerMatrix k,
    const std::vector<std::pair<int, size_t>> &c_sort,
    const std::unordered_map<int, SizeAndEntropy> &c_count,
    Rcpp::NumericVector mi,
    Rcpp::NumericVector Hmax,
    std::vector<ClusterCount> &k_counts
  )
    : k(k), c_sort(c_sort), c_count(c_count), N(k.ncol()),
      mi(mi), Hmax(Hmax), all_k_counts(k_counts), Hc(init_Hc()) {};

  void operator()(std::size_t begin, std::size_t end) {
    int c_clust;
    double Hij, Hk;
    std::vector<int> k_counts(N, 0); // number of items in each cluster of k[j]
    std::deque<int> k_counts_used; // non-zero indices in k_counts
    std::vector<int> intersects(N, 0); // number of items which occur in both k[j][l] and c[i]
    std::deque<int> intersects_used; // non-zero indices in intersects
    std::vector<size_t> cluster_sizes; // size of each cluster in k (i.e., k_counts without padding)
    std::deque<ClusterCount> my_k_counts; // multiplicity of values in cluster_sizes

    for (std::size_t j = begin; j < end; j++) {
      // generate counts for test partition k[j]
      // Rcpp::Rcerr << "generating cluster sizes for test partition " << j << std::endl;
      auto kj = k.row(j);
      k_counts_used.clear();
      for (size_t i = 0; i < N; ++i) {
        if (++k_counts[kj[i]] == 1) k_counts_used.push_back(kj[i]);
      }
      // Rcpp::Rcerr << "found " << k_counts_used.size() << " clusters" << std::endl;

      // Rcpp::Rcerr << "generating size counts" << std::endl;
      cluster_sizes.clear();
      cluster_sizes.reserve(k_counts_used.size());
      for (const auto j : k_counts_used) {
        // if (k_counts[j] == 0) Rcpp::Rcerr << "wtf, cluster of size 0??" << std::endl;
        cluster_sizes.push_back(k_counts[j]);
      }
      // Rcpp::Rcerr << "sorting size counts" << std::endl;
      std::sort(cluster_sizes.begin(), cluster_sizes.end());

      // Rcpp::Rcerr << "filling my_k_counts" << std::endl;
      // size_t pre_size = my_k_counts.size();
      Hk = 0.0;
      auto cs_i = cluster_sizes.begin();
      auto cs_end = cluster_sizes.end();
      while (cs_i != cs_end) {
        auto cs_start = cs_i;
        while (cs_i != cs_end && *cs_i == *cs_start) ++cs_i;
        size_t num = cs_i - cs_start;
        my_k_counts.push_back({*cs_start, num, j});
        Hk -= num * double(*cs_start) / N * log(double(*cs_start) / N);
      }
      // Rcpp::Rcerr << "found k clusters with " << my_k_counts.size() - pre_size
      //             << " different sizes" << std::endl;

      // Rcpp::Rcerr << "calculating MI" << std::endl;
      auto ci = c_sort.begin();
      auto c_end = c_sort.end();
      while(ci != c_end) {
        // Rcpp::Rcerr << "c cluster " << ci->first << std::endl;
        // Rcpp::Rcerr << " - finding intersects" << std::endl;
        c_clust = ci->first;
        while(ci != c_end && ci->first == c_clust) {
          if (++intersects[kj[ci->second]] == 1)
            intersects_used.push_back(kj[ci->second]);
          ++ci;
        }
        // Rcpp::Rcerr << " - sorting intersects" << std::endl;
        // std::sort(intersects.begin(), intersects.end());
        while (intersects_used.size() > 0) {
          const auto i_i = intersects_used.front();
          // Rcpp::Rcerr << "  - calculating MI for intersect with k cluster " << *i_i << std::endl;
           // Rcpp::Rcerr << "  - k cluster " << *i_i << " has " << (i_i - i_start) << " items in common with c cluster " << ci->first << std::endl;
          Hij = double(intersects[i_i]) / N;
          mi[j] += Hij * log(Hij / c_count.at(c_clust).H / double(k_counts[i_i]) * N);
          intersects_used.pop_front();
          intersects[i_i] = 0;
        }

      }
      // Rcpp::Rcerr << "calculating Hmax" << std::endl;
      Hmax[j] = Hk > Hc ? Hk : Hc;
      // Rcpp::Rcerr << "emptying k_counts" << std::endl;
      for (const auto i : k_counts_used) {
        k_counts[i] = 0;
      }
    }
    // Rcpp::Rcerr << "writing my_k_counts to all_k_counts" << std::endl;
    std::lock_guard<std::mutex> lock(mutex);
    all_k_counts.reserve(all_k_counts.size() + my_k_counts.size());
    all_k_counts.insert(all_k_counts.end(), my_k_counts.begin(), my_k_counts.end());
  }
};

class AdjustedMutualInformationWorker2 : public RcppParallel::Worker {
  RcppParallel::RVector<double> emi;
  const std::vector<ClusterCount> &k_counts;
  const std::vector<ClusterCount> c_counts;
  const std::vector<size_t> cum_calc;
  const size_t N;
  const size_t n_shard;
  const std::vector<double> lfact;
  std::mutex mutex;

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
  };

  std::vector<ClusterCount> init_c_counts(
      const std::unordered_map<int, SizeAndEntropy> &c_clust
  ) {
    std::vector<size_t> sizes;
    sizes.reserve(c_clust.size());
    for (auto c : c_clust) {
      sizes.push_back(c.second.n);
    }
    std::sort(sizes.begin(), sizes.end());

    std::vector<ClusterCount> c_counts;
    auto size_i = sizes.begin();
    auto size_end = sizes.end();
    while (size_i != size_end) {
      auto size_start = size_i;
      while (size_i != size_end && *size_i == *size_start) ++size_i;
      c_counts.push_back({*size_start, size_t(size_i - size_start), 0});
    }
    return c_counts;
  }

  std::vector<size_t> cum2(const std::vector<ClusterCount> &x,
                           const std::vector<ClusterCount> &y) {
    // Rcpp::Rcerr << "accumulating cluster count lists of length " << x.size()
    //             << " and " << y.size()
    //             << std::endl;
    std::vector<size_t> cum;
    // heuristic guess for the size; it must be between max(|x|, |y|) and |x| + |y|
    cum.reserve(std::max(x.size(), y.size()) + std::min(x.size(), y.size()) / 2);

    size_t sum = 0, sumxo = 0, sumyo = 0, sumb = 0;
    auto xi = x.begin();
    auto yi = y.begin();
    auto xend = x.end();
    auto yend = y.end();
    while (xi != xend || yi != yend) {
      if (xi == xend) {
        // Rcpp::Rcerr << "case 1a: xi == xend, y_size == " << yi->size << std::flush;
        sumyo += yi->size;
        sum += sumb + sumxo;
        size_t ysize = yi->size;
        while (yi != yend && yi->size == ysize) ++yi;
      } else if (yi == yend || xi->size < yi->size) {
        if (yi == yend) {
          // Rcpp::Rcerr << "case 2a: yi == yend, x_size == " << xi->size << std::flush;
        } else {
          // Rcpp::Rcerr << "case 2b: x_size == " << xi->size
          //             << " < y_size == " << yi->size
          //             << std::flush;
        }
        sumxo += xi->size;
        sum += sumb + sumyo;
        size_t xsize = xi->size;
        while (xi != xend && xi->size == xsize) ++xi;
      } else if (yi->size < xi->size) {
        // Rcpp::Rcerr << "case 1b: y_size == " << yi->size
        //             << " < x_size == " << xi->size
        //             << std::flush;
        sumyo += yi->size;
        sum += sumb + sumxo;
        size_t ysize = yi->size;
        while (yi != yend && yi->size == ysize) ++yi;
      } else {
        // Rcpp::Rcerr << "case 3: x_size == " << xi->size
        //             << " == y_size == " << yi->size
        //             << std::flush;
        sumb += xi->size;
        sum += sumb + sumxo + sumyo;
        size_t xsize = xi->size;
        while (xi != xend && xi->size == xsize) ++xi;
        size_t ysize = yi->size;
        while (yi != yend && yi->size == ysize) ++yi;
      }
      // Rcpp::Rcerr << "; sumb=" << sumb
      //             << " sumxo=" << sumxo
      //             << " sumyo=" << sumyo
      //             << " sum=" << sum
      //             << std::endl;
      cum.push_back(sum);
    }
    // Rcpp::Rcerr << "cumulative counts list has length " << cum.size()
    //             << std::endl;
    return cum;
  }

  double emi_term(size_t ai, size_t bi, double precalc) {
    // Rcpp::Rcerr << "emi_term for ai=" << ai
    //             << ", bi=" << bi
    //             << ", precalc=" << precalc
    //             << std::endl;
    size_t nijmin = 1;
    if (ai + bi > N) nijmin = ai + bi - N;
    double E = NEGINF;
    precalc += lfact[bi] + lfact[N - bi];
    double logN = log(N);
    double logab = log(ai) + log(bi);
    for (size_t nij = bi + 1; nij-- > nijmin;) {
      // Rcpp::Rcerr << "nij=" << nij << std::endl;
      double lognij = log(nij);
      // Rcpp::Rcerr << "log(nij)=" << lognij << std::endl;
      double det = logN - logab + lognij;
      // Rcpp::Rcerr << "det=" << det << std::endl;
      double all_but_log = lognij - logN + precalc
        - lfact[nij] - lfact[ai - nij]
        - lfact[bi - nij] - lfact[N - ai - bi + nij];
        // Rcpp::Rcerr << "all_but_log=" << all_but_log << std::endl;
      if (det > 0) {
        E = log_plus(E, all_but_log + log(det));
      } else if (det < 0) {
        E = log_minus(E, all_but_log + log(-det));
      }
      // Rcpp::Rcerr << "E=" << E << std::endl;
    }
    // Rcpp::Rcerr << ", Eij=" << E << std::flush;
    return E;
  }

public:

  AdjustedMutualInformationWorker2(
    Rcpp::NumericVector emi,
    const std::vector<ClusterCount> &k_counts,
    const std::unordered_map<int, SizeAndEntropy> &c_clust,
    const size_t N,
    const size_t n_shard
  ) :
    emi(emi),
    k_counts(k_counts),
    c_counts(init_c_counts(c_clust)),
    cum_calc(cum2(c_counts, k_counts)),
    N(N),
    n_shard(n_shard),
    lfact(init_lfact()) {};

  void operator()(size_t begin, size_t end) {
    // which calculations is this shard supposed to do?
    // Rcpp::Rcerr << "determining number of calculations" << std::endl;
    size_t n_calc = cum_calc.back();
    size_t min_calc = size_t((double)n_calc / (double)n_shard * (double)begin);
    size_t max_calc = size_t((double)n_calc / (double)n_shard * double(end));
    // Rcpp::Rcerr << "will perform calculations " << min_calc
                // << " to " << max_calc
                // << " of " << n_calc
                // << std::endl;
    std::vector<double> log_emi(emi.size(), NEGINF);

    auto calc_i = cum_calc.rbegin();
    auto calc_end = cum_calc.rend();
    auto cc_i = c_counts.rbegin();
    auto cc_end = c_counts.rend();
    auto kc_i = k_counts.rbegin();
    auto kc_end = k_counts.rend();
    while (calc_i != calc_end && kc_i != kc_end && cc_i != cc_end) {
      // size of the larger clusters; might be present in c, k, or both
      size_t asize;
      asize = kc_i->size < cc_i->size ? cc_i->size : kc_i->size;
      // Rcpp::Rcerr << "larger cluster size " << asize
      //             << " (calculations up to " << *calc_i
      //             << ")" << std::endl;
      // check if we are in the range we need to calculate
      if (*calc_i > min_calc && *calc_i <= max_calc) {
        // just do this once
        // Rcpp::Rcerr << "beginning calculations for larger cluster size " << asize
        //                << " (calculations up to " << *calc_i
        //                << ")" << std::endl;
        // << std::flush;
        double a_part = lfact[asize] + lfact[N - asize] - lfact[N];
        if (cc_i->size == kc_i->size) {
          // Rcpp::Rcerr << " (case 1: both)" << std::endl;
          auto cc_j = cc_i;
          auto kc_j = kc_i;
          do {
            size_t bsize;
            if (cc_j == cc_end) bsize = kc_j->size;
            else if (kc_j == kc_end) bsize = cc_j->size;
            else bsize = std::min(cc_j->size, kc_j->size);
            // Rcpp::Rcerr << "smaller cluster size " << bsize << std::flush;
            double Eij = emi_term(asize, bsize, a_part);
            // Rcpp::Rcerr << ", Eij=" << Eij << std::endl;
            double cEij;
            // add terms for cc_i vs kc_j
            if (kc_j != kc_end && kc_j->size == bsize) {
              // Rcpp::Rcerr << "adding terms for cc_i vs kc_j..." << std::endl;
              cEij = log(cc_i->n) + Eij;
              while (kc_j != kc_end && kc_j->size == bsize) {
                log_add_to(log_emi[kc_j->j], cEij + log(kc_j->n));
                // Rcpp::Rcerr << "log(emi[" << kc_j->j
                            // << "]) = " << log_emi[kc_j->j] << std::endl;
                ++kc_j;
              }
              // Rcpp::Rcerr << "done" << std::endl;
            }
            // add terms for kc_i vs cc_j
            if (cc_j != cc_end && cc_j->size == bsize) {
              if (asize != bsize) {
                // Rcpp::Rcerr << "adding terms for kc_i vs cc_j..." << std::flush;
                auto kc_k = kc_i;
                cEij = log(cc_j->n) + Eij;
                while (kc_k != kc_end && kc_k->size == asize) {
                  log_add_to(log_emi[kc_k->j], cEij + log(kc_k->n));
                  // Rcpp::Rcerr << "log(emi[" << kc_k->j
                  //             << "]) = " << log_emi[kc_k->j] << std::endl;
                  ++kc_k;
                }
                // Rcpp::Rcerr << "done" << std::endl;
              }
              ++cc_j;
            }
          } while (cc_j != cc_end || kc_j != kc_end);
        } else if (cc_i->size == asize) {
          // Rcpp::Rcerr << " (case 2: c only)" << std::endl;
          auto kc_j = kc_i;
          while (kc_j != kc_end && kc_j->size >= asize) ++kc_j;
          while (kc_j != kc_end) {
            size_t bsize = kc_j->size;
            // Rcpp::Rcerr << "smaller cluster size " << bsize << std::flush;
            double cEij = log(cc_i->n) + emi_term(asize, bsize, a_part);
            // Rcpp::Rcerr << ", cEij=" << cEij << std::endl;
            while(kc_j != kc_end && kc_j->size == bsize) {
              log_add_to(log_emi[kc_j->j], cEij + log(kc_j->n));
              ++kc_j;
            }
          }
        } else {
          // Rcpp::Rcerr << " (case 3: k only)" << std::endl;
          auto cc_j = cc_i;
          while (cc_j != cc_end && cc_j->size >= asize) ++cc_j;
          while (cc_j != cc_end) { //cannot be equal
            // Rcpp::Rcerr << "smaller cluster size " << cc_j->size << std::flush;
            double cEij = log(cc_j->n) + emi_term(asize, cc_j->size, a_part);
            // Rcpp::Rcerr << ", cEij=" << cEij << std::endl;
            auto kc_k = kc_i;
            while (kc_k != kc_end && kc_k->size == asize) {
              log_add_to(log_emi[kc_k->j], cEij + log(kc_k->n));
              ++kc_k;
            }
            ++cc_j;
          }
        }
      }
      if (cc_i != cc_end && cc_i->size == asize) ++cc_i;
      while (kc_i != kc_end && kc_i->size == asize) ++kc_i;
      ++calc_i;
    }
    // Rcpp::Rcerr << "remaining calculations: " << *calc_i << std::endl;
    std::lock_guard<std::mutex> lock(mutex);
    for (size_t i = 0; i < emi.size(); ++i) {
      log_add_to(emi[i], log_emi[i]);
    }
  }

};

//' @rdname mutual_information
//' @export
 // [[Rcpp::export]]
 Rcpp::DataFrame adjusted_mutual_information(
     const Rcpp::IntegerMatrix k,
     const Rcpp::IntegerVector c,
     int threads = 1L
 ) {
   size_t N = c.size(), m = k.nrow();
   if (N != (size_t)k.ncol())
     Rcpp::stop("test clusters 'k' (%d) and true clusters 'c' (%d) must have"
                  " the same number of objects.", k.ncol(), N);
   Rcpp::NumericVector mi(m, 0.0), Hmax(m, 0.0), emi(m, R_NegInf);
   std::vector<std::pair<int, size_t>> c_sort;
   std::unordered_map<int, SizeAndEntropy> c_count;
   initialize_c_counts(c, c_sort, c_count, N);

   std::vector<ClusterCount> k_counts;
   // Rcpp::Rcerr << "constructing AdjustedMutualInformationWorker1" << std::endl;
   AdjustedMutualInformationWorker1 worker1(k, c_sort, c_count, mi, Hmax, k_counts);
   // Rcpp::Rcerr << "running AdjustedMutualInformationWorker1()" << std::endl;
   if (threads == 1) {
     worker1(0, m);
   } else {
     RcppParallel::parallelFor(0, m, worker1, 1, threads);
   }
   // Rcpp::Rcerr << "finished AdjustedMutualInformationWorker1()" << std::endl;

   // Rcpp::Rcerr << "sorting k_counts" << std::endl;
   std::sort(k_counts.begin(), k_counts.end());
   // Rcpp::Rcerr << "finished sorting k_counts" << std::endl;

   // Rcpp::Rcerr << "constructing AdjustedMutualInformationWorker2" << std::endl;
   AdjustedMutualInformationWorker2 worker2(emi, k_counts, c_count, N, threads);
   // Rcpp::Rcerr << "running AdjustedMutualInformationWorker2()" << std::endl;
   if (threads == 1) {
     worker2(0, 1);
   } else {
     RcppParallel::parallelFor(0, threads, worker2, 1, threads);
   }
   // Rcpp::Rcerr << "finished AdjustedMutualInformationWorker2()" << std::endl;
   emi = exp(emi);
   Rcpp::NumericVector ami = (mi - emi) / (Hmax - emi);
   auto out = Rcpp::DataFrame::create(
     Rcpp::Named("MI") = mi,
     Rcpp::Named("EMI") = emi,
     Rcpp::Named("AMI") = ami
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
#endif //OPTIMOTU_R
