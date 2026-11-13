// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifdef OPTIMOTU_R

#include "PackedAaSequenceSet.h"
#include "AaConfig.h"
#include "EvoConfig.h"
#include "config.h"
#include <Rcpp.h>
#include <RcppParallel.h>
#include <algorithm>
#include <cstdint>
#include <mutex>
#include <unordered_map>
#include <vector>

namespace {

constexpr int kStarCutoff = 16;

uint64_t undirected_key(int a, int b) {
  const uint32_t lo = static_cast<uint32_t>(std::min(a, b));
  const uint32_t hi = static_cast<uint32_t>(std::max(a, b));
  return (static_cast<uint64_t>(hi) << 32) | lo;
}

struct AaStarWorker : public RcppParallel::Worker {
  const PackedAaSequenceSet &pss;
  const std::vector<int> &queries;
  const std::vector<std::vector<int>> &partners;
  int min_overlap;
  bool ignore_gaps;
  const float *weight_lut;
  const EvoModel *evo;
  std::unordered_map<uint64_t, double> *dists;
  std::unordered_map<uint64_t, int> *n_ok;
  std::mutex mutex;

  AaStarWorker(
    const PackedAaSequenceSet &pss,
    const std::vector<int> &queries,
    const std::vector<std::vector<int>> &partners,
    int min_overlap,
    bool ignore_gaps,
    const float *weight_lut,
    const EvoModel *evo,
    std::unordered_map<uint64_t, double> *dists,
    std::unordered_map<uint64_t, int> *n_ok
  ) : pss(pss),
      queries(queries),
      partners(partners),
      min_overlap(min_overlap),
      ignore_gaps(ignore_gaps),
      weight_lut(weight_lut),
      evo(evo),
      dists(dists),
      n_ok(n_ok) {}

  void operator()(std::size_t begin, std::size_t end) override {
    std::unordered_map<uint64_t, double> local;
    std::unordered_map<uint64_t, int> local_ok;
    for (std::size_t q = begin; q < end; ++q) {
      const int qi = queries[q];
      const auto &part = partners[static_cast<std::size_t>(qi)];
      const int n_j = static_cast<int>(part.size());
      std::vector<uint8_t> pos_major;
      pss.fill_pos_major_idx(part.data(), n_j, pos_major);
      std::vector<double> d(static_cast<std::size_t>(n_j));
      std::vector<int> n_ok(static_cast<std::size_t>(n_j));
      if (evo != nullptr) {
        pss.score_query_evo(
          qi,
          pos_major.data(),
          part.data(),
          n_j,
          n_j,
          min_overlap,
          ignore_gaps,
          *evo,
          d.data(),
          n_ok.data()
        );
      } else {
        pss.score_query_vs_pos_major(
          qi,
          pos_major.data(),
          part.data(),
          n_j,
          n_j,
          min_overlap,
          ignore_gaps,
          weight_lut,
          1.0,
          false,
          d.data(),
          n_ok.data()
        );
      }
      for (int t = 0; t < n_j; ++t) {
        const uint64_t key = undirected_key(qi, part[static_cast<std::size_t>(t)]);
        local[key] = d[static_cast<std::size_t>(t)];
        local_ok[key] = n_ok[static_cast<std::size_t>(t)];
      }
#ifdef OPTIMOTU_R
      RcppThread::checkUserInterrupt();
#endif
    }
    std::lock_guard<std::mutex> lock(mutex);
    for (const auto &kv : local) {
      (*dists)[kv.first] = kv.second;
    }
    for (const auto &kv : local_ok) {
      (*n_ok)[kv.first] = kv.second;
    }
  }
};

}  // namespace

// [[Rcpp::export]]
Rcpp::DataFrame seq_distmx_pairs_internal(
  Rcpp::CharacterVector seq,
  Rcpp::List dist_config,
  Rcpp::IntegerMatrix pairs,
  Rcpp::List parallel_config,
  int verbose = 0
) {
  (void)verbose;
  if (!dist_config.inherits("optimotu_dist_config")) {
    OPTIMOTU_STOP("'dist_config' must be of class 'optimotu_dist_config'");
  }
  const std::string method =
    element_as_string(dist_config, "method", "dist_config");
  if (method != "aa" && method != "evo") {
    OPTIMOTU_STOP(
      "seq_distmx pairs= is implemented for dist_aa() and dist_evo()"
    );
  }
  if (pairs.ncol() != 2) {
    OPTIMOTU_STOP("`pairs` must have two columns");
  }
  auto seq_str = Rcpp::as<std::vector<std::string>>(seq);
  const int n_seq = static_cast<int>(seq_str.size());
  const int n_pair = pairs.nrow();
  PackedAaSequenceSet pss(seq_str);
  std::vector<float> lut_storage;
  EvoModel evo_store;
  const EvoModel *evo = nullptr;
  const float *weight_lut = nullptr;
  if (method == "evo") {
    evo_store = evo_model_from_config(dist_config);
    evo = &evo_store;
  } else {
    weight_lut = aa_weight_lut_from_config(dist_config, lut_storage);
  }
  const int min_overlap = element_as_int(dist_config, "min_overlap", "dist_config");
  const bool ignore_gaps = element_as_bool(dist_config, "ignore_gaps", "dist_config");
  int threads = 1;
  if (parallel_config.containsElementNamed("threads")) {
    threads = element_as_int(parallel_config, "threads", "parallel_config");
  }
  if (threads < 1) {
    threads = 1;
  }

  std::vector<int> deg(static_cast<std::size_t>(n_seq), 0);
  std::vector<int> a(static_cast<std::size_t>(n_pair));
  std::vector<int> b(static_cast<std::size_t>(n_pair));
  for (int p = 0; p < n_pair; ++p) {
    int i = pairs(p, 0) - 1;
    int j = pairs(p, 1) - 1;
    if (i < 0 || j < 0 || i >= n_seq || j >= n_seq) {
      OPTIMOTU_STOP("`pairs` contains an index outside 1:length(seq)");
    }
    if (i == j) {
      OPTIMOTU_STOP("`pairs` cannot contain self-comparisons");
    }
    a[static_cast<std::size_t>(p)] = i;
    b[static_cast<std::size_t>(p)] = j;
    ++deg[static_cast<std::size_t>(i)];
    ++deg[static_cast<std::size_t>(j)];
  }

  std::vector<std::vector<int>> partners(static_cast<std::size_t>(n_seq));
  std::vector<char> seen_star(static_cast<std::size_t>(n_seq), 0);
  std::unordered_map<uint64_t, double> dist_map;
  std::unordered_map<uint64_t, int> n_ok_map;
  std::vector<int> leftover_i;
  std::vector<int> leftover_j;
  leftover_i.reserve(static_cast<std::size_t>(n_pair));
  leftover_j.reserve(static_cast<std::size_t>(n_pair));

  for (int p = 0; p < n_pair; ++p) {
    const int i = a[static_cast<std::size_t>(p)];
    const int j = b[static_cast<std::size_t>(p)];
    const uint64_t key = undirected_key(i, j);
    if (dist_map.find(key) != dist_map.end()) {
      continue;
    }
    dist_map[key] = -1.0;
    const int di = deg[static_cast<std::size_t>(i)];
    const int dj = deg[static_cast<std::size_t>(j)];
    int q = i;
    int t = j;
    if (dj > di || (dj == di && j < i)) {
      q = j;
      t = i;
    }
    if (deg[static_cast<std::size_t>(q)] >= kStarCutoff) {
      partners[static_cast<std::size_t>(q)].push_back(t);
      seen_star[static_cast<std::size_t>(q)] = 1;
    } else {
      leftover_i.push_back(i);
      leftover_j.push_back(j);
    }
  }

  std::vector<int> queries;
  queries.reserve(static_cast<std::size_t>(n_seq));
  for (int i = 0; i < n_seq; ++i) {
    if (seen_star[static_cast<std::size_t>(i)]) {
      queries.push_back(i);
    }
  }

  if (!queries.empty()) {
    AaStarWorker worker(
      pss, queries, partners, min_overlap, ignore_gaps, weight_lut, evo,
      &dist_map, &n_ok_map
    );
    const int n_threads = std::min(threads, static_cast<int>(queries.size()));
    if (n_threads <= 1) {
      worker(0, queries.size());
    } else {
      RcppParallel::parallelFor(0, queries.size(), worker, 1, n_threads);
    }
  }

  for (std::size_t k = 0; k < leftover_i.size(); ++k) {
    const int i = leftover_i[k];
    const int j = leftover_j[k];
    const auto scored = evo == nullptr
      ? pss.success_score_and_dist(i, j, min_overlap, ignore_gaps, weight_lut)
      : pss.score_evo_pair(i, j, min_overlap, ignore_gaps, *evo);
    const uint64_t key = undirected_key(i, j);
    dist_map[key] = std::get<2>(scored);
    n_ok_map[key] = std::get<1>(scored);
  }

  Rcpp::IntegerVector seq1_out(n_pair);
  Rcpp::IntegerVector seq2_out(n_pair);
  Rcpp::NumericVector score1_out(n_pair, 0.0);
  Rcpp::NumericVector score2_out(n_pair, 0.0);
  Rcpp::NumericVector dist1_out(n_pair, 0.0);
  Rcpp::NumericVector dist2_out(n_pair);
  for (int p = 0; p < n_pair; ++p) {
    const int i = a[static_cast<std::size_t>(p)];
    const int j = b[static_cast<std::size_t>(p)];
    seq1_out[p] = i + 1;
    seq2_out[p] = j + 1;
    const uint64_t key = undirected_key(i, j);
    dist2_out[p] = dist_map[key];
    score2_out[p] = n_ok_map[key];
  }
  return Rcpp::DataFrame::create(
    Rcpp::Named("seq_idx1") = seq1_out,
    Rcpp::Named("seq_idx2") = seq2_out,
    Rcpp::Named("score1") = score1_out,
    Rcpp::Named("score2") = score2_out,
    Rcpp::Named("dist1") = dist1_out,
    Rcpp::Named("dist2") = dist2_out
  );
}

#endif
