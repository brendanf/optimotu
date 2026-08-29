// SPDX-FileCopyrightText: 2025 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifdef OPTIMOTU_R

#include "cluster_measures.h"
#include "optimotu.h"
#include <RcppParallel.h>
#include <cmath>
#include <deque>
#include <limits>
#include <mutex>
#include <unordered_map>
#include <algorithm>

#define NEGINF (-std::numeric_limits<double>::infinity())

namespace {

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

struct ClusterCount {
  size_t size;
  size_t n;
  size_t j;
  bool operator<(const ClusterCount &x) const { return size < x.size; }
};

void initialize_c_counts(
  const Rcpp::IntegerVector &c,
  std::vector<std::pair<int, size_t>> &c_sort,
  std::unordered_map<int, size_t> &c_count,
  size_t N
) {
  c_sort.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    c_sort.emplace_back(c[i], i);
    auto it = c_count.find(c[i]);
    if (it == c_count.end()) {
      c_count.emplace(c[i], 1);
    } else {
      ++it->second;
    }
  }
  std::sort(c_sort.begin(), c_sort.end());
}

void initialize_c_entropy(
  const Rcpp::IntegerVector &c,
  std::vector<std::pair<int, size_t>> &c_sort,
  std::unordered_map<int, SizeAndEntropy> &c_count,
  size_t N
) {
  double Ninv = 1.0 / N;
  c_sort.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    c_sort.emplace_back(c[i], i);
    auto it = c_count.find(c[i]);
    if (it == c_count.end()) {
      c_count.emplace(c[i], SizeAndEntropy{Ninv, 1});
    } else {
      ++it->second.n;
      it->second.H = it->second.n * Ninv;
    }
  }
  std::sort(c_sort.begin(), c_sort.end());
}

void confusion_row(
  const int *kj,
  size_t N,
  const std::vector<std::pair<int, size_t>> &c_sort,
  const std::unordered_map<int, size_t> &c_count,
  double &tp,
  double &fp,
  double &fn,
  double &tn
) {
  const size_t Npairs = N * (N - 1) / 2;
  std::unordered_map<int, size_t> k_count;
  std::unordered_map<int, size_t> intersects;
  for (size_t i = 0; i < N; ++i) {
    auto it = k_count.find(kj[i]);
    if (it == k_count.end()) {
      k_count.emplace(kj[i], 1);
    } else {
      ++it->second;
    }
  }
  size_t tpj = 0, fpj = 0, fnj = 0;
  auto ci = c_sort.begin();
  auto c_end = c_sort.end();
  while (ci != c_end) {
    int c_clust = ci->first;
    intersects.clear();
    while (ci != c_end && ci->first == c_clust) {
      int k_clust = kj[ci->second];
      auto it = intersects.find(k_clust);
      if (it == intersects.end()) {
        intersects.emplace(k_clust, 1);
      } else {
        ++it->second;
      }
      ++ci;
    }
    for (const auto &counti : intersects) {
      tpj += counti.second * (counti.second - 1);
      fpj += counti.second * (k_count.at(counti.first) - counti.second);
      fnj += counti.second * (c_count.at(c_clust) - counti.second);
    }
  }
  tp = tpj / 2.0;
  fp = fpj / 2.0;
  fn = fnj / 2.0;
  tn = static_cast<double>(Npairs) - tp - fp - fn;
}

double fmeasure_row(
  const int *kj,
  size_t N,
  const std::vector<std::pair<int, size_t>> &c_sort,
  const std::unordered_map<int, size_t> &c_count
) {
  std::unordered_map<int, size_t> k_count;
  std::unordered_map<int, size_t> intersects;
  for (size_t i = 0; i < N; ++i) {
    auto it = k_count.find(kj[i]);
    if (it == k_count.end()) {
      k_count.emplace(kj[i], 1);
    } else {
      ++it->second;
    }
  }
  double result = 0.0;
  auto ci = c_sort.begin();
  auto c_end = c_sort.end();
  while (ci != c_end) {
    int c_clust = ci->first;
    intersects.clear();
    while (ci != c_end && ci->first == c_clust) {
      int k_clust = kj[ci->second];
      auto it = intersects.find(k_clust);
      if (it == intersects.end()) {
        intersects.emplace(k_clust, 1);
      } else {
        ++it->second;
      }
      ++ci;
    }
    double m = 0;
    for (const auto &counti : intersects) {
      double mm = counti.second /
        static_cast<double>(c_count.at(c_clust) + k_count.at(counti.first));
      if (mm > m) m = mm;
    }
    result += c_count.at(c_clust) * m;
  }
  return result * 2.0 / N;
}

void ami_row(
  const int *kj,
  size_t N,
  const std::vector<std::pair<int, size_t>> &c_sort,
  const std::unordered_map<int, SizeAndEntropy> &c_count,
  double Hc,
  double &mi,
  double &Hmax,
  std::deque<ClusterCount> &my_k_counts,
  size_t j
) {
  std::vector<int> k_counts(N, 0);
  std::deque<int> k_counts_used;
  std::vector<int> intersects(N, 0);
  std::deque<int> intersects_used;
  for (size_t i = 0; i < N; ++i) {
    int id = kj[i];
    if (id < 0 || static_cast<size_t>(id) >= N) {
      OPTIMOTU_STOP("cluster id out of range in AMI calculation");
    }
    if (++k_counts[id] == 1) k_counts_used.push_back(id);
  }
  std::vector<size_t> cluster_sizes;
  cluster_sizes.reserve(k_counts_used.size());
  for (int id : k_counts_used) {
    cluster_sizes.push_back(k_counts[id]);
  }
  std::sort(cluster_sizes.begin(), cluster_sizes.end());
  double Hk = 0.0;
  auto cs_i = cluster_sizes.begin();
  auto cs_end = cluster_sizes.end();
  while (cs_i != cs_end) {
    auto cs_start = cs_i;
    while (cs_i != cs_end && *cs_i == *cs_start) ++cs_i;
    size_t num = cs_i - cs_start;
    my_k_counts.push_back({*cs_start, num, j});
    Hk -= num * double(*cs_start) / N * std::log(double(*cs_start) / N);
  }
  mi = 0.0;
  auto ci = c_sort.begin();
  auto c_end = c_sort.end();
  while (ci != c_end) {
    int c_clust = ci->first;
    while (ci != c_end && ci->first == c_clust) {
      int id = kj[ci->second];
      if (++intersects[id] == 1) intersects_used.push_back(id);
      ++ci;
    }
    while (!intersects_used.empty()) {
      const auto i_i = intersects_used.front();
      double Hij = double(intersects[i_i]) / N;
      mi += Hij * std::log(
        Hij / c_count.at(c_clust).H / double(k_counts[i_i]) * N
      );
      intersects_used.pop_front();
      intersects[i_i] = 0;
    }
  }
  Hmax = Hk > Hc ? Hk : Hc;
}

double mcc_from(double tp, double fp, double fn, double tn) {
  return (tp * tn - fp * fn) /
    std::sqrt(tp + fp) / std::sqrt(tp + fn) /
    std::sqrt(fn + tn) / std::sqrt(fp + tn);
}

double ri_from(double tp, double fp, double fn, double tn) {
  return (tp + tn) / (tp + fp + fn + tn);
}

double ari_from(double tp, double fp, double fn, double tn) {
  return 2.0 * (tp * tn - fp * fn) /
    ((tp + fp) * (fp + tn) + (tp + fn) * (fn + tn));
}

double fmi_from(double tp, double fp, double fn) {
  return tp / std::sqrt(tp + fp) / std::sqrt(tp + fn);
}

// 0-based index matching R strict_median(which(values == max)).
std::size_t pick_best_index(const std::vector<double> &values) {
  double best = NEGINF;
  bool any = false;
  for (double v : values) {
    if (std::isfinite(v) && (!any || v > best)) {
      best = v;
      any = true;
    }
  }
  if (!any) return values.size();
  std::vector<std::size_t> hits;
  for (std::size_t i = 0; i < values.size(); ++i) {
    if (values[i] == best) hits.push_back(i);
  }
  return hits[(hits.size() - 1) / 2];
}

class AmiWorker2 : public RcppParallel::Worker {
  std::vector<double> &emi;
  const std::vector<ClusterCount> &k_counts;
  const std::vector<ClusterCount> c_counts;
  const std::vector<size_t> cum_calc;
  const size_t N;
  const size_t n_shard;
  const std::vector<double> lfact;
  std::mutex mutex;

  std::vector<double> init_lfact() const {
    std::vector<double> v;
    v.reserve(N + 1);
    double s = 0.0;
    v.push_back(s);
    for (size_t i = 1; i <= N; ++i) {
      s += std::log(i);
      v.push_back(s);
    }
    return v;
  }

  static std::vector<ClusterCount> init_c_counts(
    const std::unordered_map<int, SizeAndEntropy> &c_clust
  ) {
    std::vector<size_t> sizes;
    sizes.reserve(c_clust.size());
    for (auto c : c_clust) sizes.push_back(c.second.n);
    std::sort(sizes.begin(), sizes.end());
    std::vector<ClusterCount> out;
    auto size_i = sizes.begin();
    auto size_end = sizes.end();
    while (size_i != size_end) {
      auto size_start = size_i;
      while (size_i != size_end && *size_i == *size_start) ++size_i;
      out.push_back({*size_start, size_t(size_i - size_start), 0});
    }
    return out;
  }

  static std::vector<size_t> cum2(
    const std::vector<ClusterCount> &x,
    const std::vector<ClusterCount> &y
  ) {
    std::vector<size_t> cum;
    cum.reserve(std::max(x.size(), y.size()) + std::min(x.size(), y.size()) / 2);
    size_t sum = 0, sumxo = 0, sumyo = 0, sumb = 0;
    auto xi = x.begin();
    auto yi = y.begin();
    auto xend = x.end();
    auto yend = y.end();
    while (xi != xend || yi != yend) {
      if (xi == xend) {
        sumyo += yi->size;
        sum += sumb + sumxo;
        size_t ysize = yi->size;
        while (yi != yend && yi->size == ysize) ++yi;
      } else if (yi == yend || xi->size < yi->size) {
        sumxo += xi->size;
        sum += sumb + sumyo;
        size_t xsize = xi->size;
        while (xi != xend && xi->size == xsize) ++xi;
      } else if (yi->size < xi->size) {
        sumyo += yi->size;
        sum += sumb + sumxo;
        size_t ysize = yi->size;
        while (yi != yend && yi->size == ysize) ++yi;
      } else {
        sumb += xi->size;
        sum += sumb + sumxo + sumyo;
        size_t xsize = xi->size;
        while (xi != xend && xi->size == xsize) ++xi;
        size_t ysize = yi->size;
        while (yi != yend && yi->size == ysize) ++yi;
      }
      cum.push_back(sum);
    }
    return cum;
  }

  double emi_term(size_t ai, size_t bi, double precalc) const {
    size_t nijmin = 1;
    if (ai + bi > N) nijmin = ai + bi - N;
    double E = NEGINF;
    precalc += lfact[bi] + lfact[N - bi];
    double logN = std::log(N);
    double logab = std::log(ai) + std::log(bi);
    for (size_t nij = bi + 1; nij-- > nijmin;) {
      double lognij = std::log(nij);
      double det = logN - logab + lognij;
      double all_but_log = lognij - logN + precalc
        - lfact[nij] - lfact[ai - nij]
        - lfact[bi - nij] - lfact[N - ai - bi + nij];
      if (det > 0) {
        E = log_plus(E, all_but_log + std::log(det));
      } else if (det < 0) {
        E = log_minus(E, all_but_log + std::log(-det));
      }
    }
    return E;
  }

public:
  AmiWorker2(
    std::vector<double> &emi,
    const std::vector<ClusterCount> &k_counts,
    const std::unordered_map<int, SizeAndEntropy> &c_clust,
    size_t N,
    size_t n_shard
  ) :
    emi(emi),
    k_counts(k_counts),
    c_counts(init_c_counts(c_clust)),
    cum_calc(cum2(c_counts, k_counts)),
    N(N),
    n_shard(n_shard),
    lfact(init_lfact()) {}

  void operator()(size_t begin, size_t end) {
    size_t n_calc = cum_calc.empty() ? 0 : cum_calc.back();
    size_t min_calc = size_t(
      (double)n_calc / (double)n_shard * (double)begin
    );
    size_t max_calc = size_t(
      (double)n_calc / (double)n_shard * double(end)
    );
    std::vector<double> log_emi(emi.size(), NEGINF);
    auto calc_i = cum_calc.rbegin();
    auto calc_end = cum_calc.rend();
    auto cc_i = c_counts.rbegin();
    auto cc_end = c_counts.rend();
    auto kc_i = k_counts.rbegin();
    auto kc_end = k_counts.rend();
    while (calc_i != calc_end && kc_i != kc_end && cc_i != cc_end) {
      size_t asize = kc_i->size < cc_i->size ? cc_i->size : kc_i->size;
      if (*calc_i > min_calc && *calc_i <= max_calc) {
        double a_part = lfact[asize] + lfact[N - asize] - lfact[N];
        if (cc_i->size == kc_i->size) {
          auto cc_j = cc_i;
          auto kc_j = kc_i;
          do {
            size_t bsize;
            if (cc_j == cc_end) bsize = kc_j->size;
            else if (kc_j == kc_end) bsize = cc_j->size;
            else bsize = std::min(cc_j->size, kc_j->size);
            double Eij = emi_term(asize, bsize, a_part);
            if (kc_j != kc_end && kc_j->size == bsize) {
              double cEij = std::log(cc_i->n) + Eij;
              while (kc_j != kc_end && kc_j->size == bsize) {
                log_add_to(log_emi[kc_j->j], cEij + std::log(kc_j->n));
                ++kc_j;
              }
            }
            if (cc_j != cc_end && cc_j->size == bsize) {
              if (asize != bsize) {
                auto kc_k = kc_i;
                double cEij = std::log(cc_j->n) + Eij;
                while (kc_k != kc_end && kc_k->size == asize) {
                  log_add_to(log_emi[kc_k->j], cEij + std::log(kc_k->n));
                  ++kc_k;
                }
              }
              ++cc_j;
            }
          } while (cc_j != cc_end || kc_j != kc_end);
        } else if (cc_i->size == asize) {
          auto kc_j = kc_i;
          while (kc_j != kc_end && kc_j->size >= asize) ++kc_j;
          while (kc_j != kc_end) {
            size_t bsize = kc_j->size;
            double cEij = std::log(cc_i->n) + emi_term(asize, bsize, a_part);
            while (kc_j != kc_end && kc_j->size == bsize) {
              log_add_to(log_emi[kc_j->j], cEij + std::log(kc_j->n));
              ++kc_j;
            }
          }
        } else {
          auto cc_j = cc_i;
          while (cc_j != cc_end && cc_j->size >= asize) ++cc_j;
          while (cc_j != cc_end) {
            double cEij = std::log(cc_j->n) + emi_term(asize, cc_j->size, a_part);
            auto kc_k = kc_i;
            while (kc_k != kc_end && kc_k->size == asize) {
              log_add_to(log_emi[kc_k->j], cEij + std::log(kc_k->n));
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
    std::lock_guard<std::mutex> lock(mutex);
    for (size_t i = 0; i < emi.size(); ++i) {
      log_add_to(emi[i], log_emi[i]);
    }
  }
};

struct ScoreWorker : public RcppParallel::Worker {
  MultipleClusterAlgorithm &algo;
  const std::size_t subset;
  const std::size_t N;
  const std::vector<std::pair<int, size_t>> &c_sort;
  const std::unordered_map<int, size_t> &c_count;
  const std::unordered_map<int, SizeAndEntropy> *c_entropy;
  const double Hc;
  const bool need_confusion;
  const bool need_ami;
  const bool need_fm;
  std::vector<double> &tp, &fp, &fn, &tn;
  std::vector<double> &mi, &Hmax;
  std::vector<double> &fm;
  std::vector<ClusterCount> &all_k_counts;
  std::mutex mutex;

  ScoreWorker(
    MultipleClusterAlgorithm &algo,
    std::size_t subset,
    std::size_t N,
    const std::vector<std::pair<int, size_t>> &c_sort,
    const std::unordered_map<int, size_t> &c_count,
    const std::unordered_map<int, SizeAndEntropy> *c_entropy,
    double Hc,
    bool need_confusion,
    bool need_ami,
    bool need_fm,
    std::vector<double> &tp,
    std::vector<double> &fp,
    std::vector<double> &fn,
    std::vector<double> &tn,
    std::vector<double> &mi,
    std::vector<double> &Hmax,
    std::vector<double> &fm,
    std::vector<ClusterCount> &all_k_counts
  ) :
    algo(algo), subset(subset), N(N), c_sort(c_sort), c_count(c_count),
    c_entropy(c_entropy), Hc(Hc),
    need_confusion(need_confusion), need_ami(need_ami), need_fm(need_fm),
    tp(tp), fp(fp), fn(fn), tn(tn), mi(mi), Hmax(Hmax), fm(fm),
    all_k_counts(all_k_counts) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<int> buf(N);
    std::deque<ClusterCount> my_k_counts;
    for (std::size_t j = begin; j < end; ++j) {
      algo.write_threshold_row(subset, static_cast<d_t>(j), buf.data());
      if (need_confusion) {
        confusion_row(
          buf.data(), N, c_sort, c_count, tp[j], fp[j], fn[j], tn[j]
        );
      }
      if (need_fm) {
        fm[j] = fmeasure_row(buf.data(), N, c_sort, c_count);
      }
      if (need_ami) {
        ami_row(
          buf.data(), N, c_sort, *c_entropy, Hc,
          mi[j], Hmax[j], my_k_counts, j
        );
      }
    }
    if (need_ami && !my_k_counts.empty()) {
      std::lock_guard<std::mutex> lock(mutex);
      all_k_counts.insert(
        all_k_counts.end(), my_k_counts.begin(), my_k_counts.end()
      );
    }
  }
};

std::vector<double> expand_scores(
  const std::vector<double> &bin_scores,
  const Rcpp::IntegerVector &threshold_order
) {
  if (threshold_order.size() == 0) return bin_scores;
  std::vector<double> out(threshold_order.size());
  for (R_xlen_t i = 0; i < threshold_order.size(); ++i) {
    int idx = threshold_order[i] - 1;
    if (idx < 0 || static_cast<size_t>(idx) >= bin_scores.size()) {
      OPTIMOTU_STOP("threshold_order index out of range");
    }
    out[i] = bin_scores[idx];
  }
  return out;
}

} // namespace

Rcpp::DataFrame find_best_threshold_from_algo(
  MultipleClusterAlgorithm &algo,
  std::size_t subset,
  const Rcpp::IntegerVector &c,
  const std::vector<std::string> &measures,
  const Rcpp::NumericVector &thresholds,
  const Rcpp::IntegerVector &threshold_order,
  int threads
) {
  const size_t N = static_cast<size_t>(c.size());
  const d_t m = algo.m;
  if (N != static_cast<size_t>(algo.subset_n(subset))) {
    OPTIMOTU_STOP(
      "true partition length does not match subset size"
    );
  }
  bool need_mcc = false, need_ri = false, need_ari = false, need_fmi = false;
  bool need_mi = false, need_ami = false, need_fm = false;
  std::vector<std::string> sorted_measures;
  for (const auto &meas : measures) {
    if (meas == "MCC") {
      need_mcc = true;
      sorted_measures.push_back(meas);
    } else if (meas == "RI") {
      need_ri = true;
      sorted_measures.push_back(meas);
    } else if (meas == "ARI") {
      need_ari = true;
      sorted_measures.push_back(meas);
    } else if (meas == "FMI") {
      need_fmi = true;
      sorted_measures.push_back(meas);
    } else if (meas == "MI") {
      need_mi = true;
    } else if (meas == "AMI") {
      need_ami = true;
    } else if (meas == "FM") {
      need_fm = true;
    } else {
      OPTIMOTU_STOP("unknown clustering measure");
    }
  }
  if (need_mi) sorted_measures.push_back("MI");
  if (need_ami) sorted_measures.push_back("AMI");
  if (need_fm) sorted_measures.push_back("FM");

  const bool need_confusion = need_mcc || need_ri || need_ari || need_fmi;
  const bool need_ami_any = need_mi || need_ami;

  std::vector<std::pair<int, size_t>> c_sort;
  std::unordered_map<int, size_t> c_count;
  initialize_c_counts(c, c_sort, c_count, N);

  std::unordered_map<int, SizeAndEntropy> c_entropy;
  std::vector<std::pair<int, size_t>> c_sort_ami;
  double Hc = 0.0;
  if (need_ami_any) {
    initialize_c_entropy(c, c_sort_ami, c_entropy, N);
    for (const auto &ci : c_entropy) {
      Hc -= ci.second.H * std::log(ci.second.H);
    }
  }

  std::vector<double> tp(m, 0), fp(m, 0), fn(m, 0), tn(m, 0);
  std::vector<double> mi(m, 0), Hmax(m, 0), fm(m, 0);
  std::vector<ClusterCount> all_k_counts;

  ScoreWorker worker(
    algo, subset, N,
    need_ami_any ? c_sort_ami : c_sort,
    c_count,
    need_ami_any ? &c_entropy : nullptr,
    Hc,
    need_confusion, need_ami_any, need_fm,
    tp, fp, fn, tn, mi, Hmax, fm, all_k_counts
  );
  if (threads <= 1) {
    worker(0, static_cast<std::size_t>(m));
  } else {
    RcppParallel::parallelFor(0, static_cast<std::size_t>(m), worker, 1, threads);
  }

  std::vector<double> ami_vals(m, NA_REAL);
  if (need_ami_any) {
    std::sort(all_k_counts.begin(), all_k_counts.end());
    std::vector<double> emi(m, NEGINF);
    AmiWorker2 worker2(emi, all_k_counts, c_entropy, N, threads <= 1 ? 1 : threads);
    if (threads <= 1) {
      worker2(0, 1);
    } else {
      RcppParallel::parallelFor(0, static_cast<std::size_t>(threads), worker2, 1, threads);
    }
    for (d_t j = 0; j < m; ++j) {
      emi[j] = std::exp(emi[j]);
      ami_vals[j] = (mi[j] - emi[j]) / (Hmax[j] - emi[j]);
    }
  }

  auto bin_to_out = [&](const std::vector<double> &bin) {
    return expand_scores(bin, threshold_order);
  };

  Rcpp::CharacterVector out_measure(sorted_measures.size());
  Rcpp::NumericVector out_threshold(sorted_measures.size(), NA_REAL);
  Rcpp::NumericVector out_value(sorted_measures.size(), NA_REAL);

  auto set_best = [&](int i, const std::vector<double> &expanded) {
    std::size_t idx = pick_best_index(expanded);
    if (idx >= expanded.size()) {
      // Match R max(..., na.rm = TRUE) on an all-NA vector.
      out_value[i] = NEGINF;
      return;
    }
    out_threshold[i] = thresholds[idx];
    out_value[i] = expanded[idx];
  };

  int i = 0;
  for (const auto &meas : sorted_measures) {
    out_measure[i] = meas;
    std::vector<double> bin(m);
    if (meas == "MCC") {
      for (d_t j = 0; j < m; ++j) bin[j] = mcc_from(tp[j], fp[j], fn[j], tn[j]);
    } else if (meas == "RI") {
      for (d_t j = 0; j < m; ++j) bin[j] = ri_from(tp[j], fp[j], fn[j], tn[j]);
    } else if (meas == "ARI") {
      for (d_t j = 0; j < m; ++j) bin[j] = ari_from(tp[j], fp[j], fn[j], tn[j]);
    } else if (meas == "FMI") {
      for (d_t j = 0; j < m; ++j) bin[j] = fmi_from(tp[j], fp[j], fn[j]);
    } else if (meas == "MI") {
      bin = mi;
    } else if (meas == "AMI") {
      bin = ami_vals;
    } else if (meas == "FM") {
      bin = fm;
    }
    set_best(i, bin_to_out(bin));
    ++i;
  }

  return Rcpp::DataFrame::create(
    Rcpp::Named("measure") = out_measure,
    Rcpp::Named("threshold") = out_threshold,
    Rcpp::Named("value") = out_value
  );
}

Rcpp::IntegerMatrix subset_matrix_via_rows(
  MultipleClusterAlgorithm &algo,
  std::size_t subset
) {
  const j_t n = algo.subset_n(subset);
  const d_t m = algo.m;
  Rcpp::IntegerMatrix out(m, n);
  std::vector<int> row(n);
  for (d_t t = 0; t < m; ++t) {
    algo.write_threshold_row(subset, t, row.data());
    for (j_t j = 0; j < n; ++j) {
      out(t, j) = row[j];
    }
  }
  return out;
}

#endif // OPTIMOTU_R
