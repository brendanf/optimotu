// SPDX-FileCopyrightText: 2026 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "profile_stats.h"
#include "optimotu.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <mutex>
#include <sstream>
#include <vector>

namespace {

struct MergeSpan {
  double start_s = 0;
  double end_s = 0;
};

struct TileTls {
  bool active = false;
  std::size_t index = 0;
  std::uint64_t pairs = 0;
  std::chrono::steady_clock::time_point start{};
  std::chrono::steady_clock::time_point pair_done{};
  bool have_pair_done = false;
  double finish_s = 0;
};

struct ProfileState {
  std::mutex mu;
  std::vector<ProfileTile> tiles;
  double worker_s = 0;
  double finalize_s = 0;
  std::chrono::steady_clock::time_point origin{};
  bool have_origin = false;
  std::atomic<double> last_tile_work_s{0};
  std::atomic<int> n_tile_marks{0};
  std::vector<MergeSpan> merges;
  double merge_complete_s = 0;
  bool have_merge_complete = false;
};

double now_s(const ProfileState &s) {
  return std::chrono::duration<double>(
             std::chrono::steady_clock::now() - s.origin)
      .count();
}

double &lock_wait_tls() {
  thread_local double seconds = 0;
  return seconds;
}

std::uint64_t &which_calls_tls() {
  thread_local std::uint64_t calls = 0;
  return calls;
}

std::uint64_t &which_sum_tls() {
  thread_local std::uint64_t sum = 0;
  return sum;
}

std::uint64_t &applied_tls() {
  thread_local std::uint64_t n = 0;
  return n;
}

std::uint64_t &ksw2_aligned_tls() {
  thread_local std::uint64_t n = 0;
  return n;
}

std::uint64_t &ksw2_skipped_tls() {
  thread_local std::uint64_t n = 0;
  return n;
}

TileTls &tile_tls() {
  thread_local TileTls t;
  return t;
}

ProfileState &state() {
  static ProfileState s;
  return s;
}

double percentile_sorted(const std::vector<double> &x, double p) {
  if (x.empty()) {
    return 0;
  }
  const double idx = p * static_cast<double>(x.size() - 1);
  const std::size_t i = static_cast<std::size_t>(std::floor(idx));
  const std::size_t j = std::min(i + 1, x.size() - 1);
  const double f = idx - static_cast<double>(i);
  return x[i] * (1.0 - f) + x[j] * f;
}

void write_tile_tsv(const std::vector<ProfileTile> &tiles) {
  const char *path = std::getenv("OPTIMOTU_PROFILE_TILES");
  if (path == nullptr || path[0] == '\0') {
    return;
  }
  std::ifstream check(path);
  const bool empty = !check.good() ||
                     check.peek() == std::ifstream::traits_type::eof();
  check.close();
  std::ofstream out(path, std::ios::app);
  if (!out) {
    return;
  }
  if (empty) {
    out << "tile_index\tpairs\twall_s\tpair_s\tfinish_s\twhich_mean\t"
           "applied\tksw2_aligned\tksw2_skipped\n";
  }
  for (const auto &t : tiles) {
    const double wall = t.setup_s + t.pair_s + t.finish_s;
    const double which_mean = t.which_calls > 0
        ? static_cast<double>(t.which_sum) /
              static_cast<double>(t.which_calls)
        : 0;
    out << t.tile_index << '\t' << t.pairs << '\t' << wall << '\t'
        << t.pair_s << '\t' << t.finish_s << '\t' << which_mean << '\t'
        << t.applied << '\t' << t.ksw2_aligned << '\t' << t.ksw2_skipped
        << '\n';
  }
}

} // namespace

bool optimotu_profile_enabled() {
  static const bool on = []() {
    const char *e = std::getenv("OPTIMOTU_PROFILE");
    return e != nullptr && e[0] != '\0' && e[0] != '0';
  }();
  return on;
}

void optimotu_profile_reset() {
  if (!optimotu_profile_enabled()) {
    return;
  }
  auto &s = state();
  std::lock_guard<std::mutex> lock(s.mu);
  s.tiles.clear();
  s.worker_s = 0;
  s.finalize_s = 0;
  s.have_origin = false;
  s.last_tile_work_s.store(0);
  s.n_tile_marks.store(0);
  s.merges.clear();
  s.merge_complete_s = 0;
  s.have_merge_complete = false;
}

void optimotu_profile_add_tile(const ProfileTile &tile) {
  auto &s = state();
  std::lock_guard<std::mutex> lock(s.mu);
  s.tiles.push_back(tile);
}

void optimotu_profile_add_lock_wait(double seconds) {
  lock_wait_tls() += seconds;
}

double optimotu_profile_take_lock_wait() {
  double seconds = lock_wait_tls();
  lock_wait_tls() = 0;
  return seconds;
}

void optimotu_profile_add_which(std::uint32_t n_subsets) {
  ++which_calls_tls();
  which_sum_tls() += n_subsets;
}

void optimotu_profile_take_which(
    std::uint64_t &calls,
    std::uint64_t &sum) {
  calls = which_calls_tls();
  sum = which_sum_tls();
  which_calls_tls() = 0;
  which_sum_tls() = 0;
}

void optimotu_profile_add_applied(std::uint64_t n) {
  applied_tls() += n;
}

void optimotu_profile_take_applied(std::uint64_t &n) {
  n = applied_tls();
  applied_tls() = 0;
}

void optimotu_profile_add_ksw2(std::uint64_t aligned, std::uint64_t skipped) {
  ksw2_aligned_tls() += aligned;
  ksw2_skipped_tls() += skipped;
}

void optimotu_profile_take_ksw2(
    std::uint64_t &aligned,
    std::uint64_t &skipped) {
  aligned = ksw2_aligned_tls();
  skipped = ksw2_skipped_tls();
  ksw2_aligned_tls() = 0;
  ksw2_skipped_tls() = 0;
}

void optimotu_profile_set_worker(double seconds) {
  auto &s = state();
  std::lock_guard<std::mutex> lock(s.mu);
  s.worker_s = seconds;
}

void optimotu_profile_set_finalize(double seconds) {
  auto &s = state();
  std::lock_guard<std::mutex> lock(s.mu);
  s.finalize_s = seconds;
}

void optimotu_profile_mark_origin() {
  if (!optimotu_profile_enabled()) {
    return;
  }
  auto &s = state();
  std::lock_guard<std::mutex> lock(s.mu);
  s.origin = std::chrono::steady_clock::now();
  s.have_origin = true;
  s.last_tile_work_s.store(0);
  s.n_tile_marks.store(0);
  s.merges.clear();
  s.merge_complete_s = 0;
  s.have_merge_complete = false;
}

double optimotu_profile_now_s() {
  if (!optimotu_profile_enabled()) {
    return 0;
  }
  auto &s = state();
  if (!s.have_origin) {
    return 0;
  }
  return now_s(s);
}

void optimotu_profile_mark_tile_work_done() {
  if (!optimotu_profile_enabled()) {
    return;
  }
  auto &s = state();
  if (!s.have_origin) {
    return;
  }
  const double t = now_s(s);
  double prev = s.last_tile_work_s.load(std::memory_order_relaxed);
  while (t > prev &&
         !s.last_tile_work_s.compare_exchange_weak(
             prev, t, std::memory_order_relaxed)) {
  }
  s.n_tile_marks.fetch_add(1, std::memory_order_relaxed);
}

void optimotu_profile_add_merge_span(double start_s, double end_s) {
  if (!optimotu_profile_enabled()) {
    return;
  }
  auto &s = state();
  std::lock_guard<std::mutex> lock(s.mu);
  s.merges.push_back({start_s, end_s});
}

void optimotu_profile_mark_merge_complete() {
  if (!optimotu_profile_enabled()) {
    return;
  }
  auto &s = state();
  std::lock_guard<std::mutex> lock(s.mu);
  if (!s.have_origin) {
    return;
  }
  s.merge_complete_s = now_s(s);
  s.have_merge_complete = true;
}

void optimotu_profile_tile_enter(std::size_t index, std::uint64_t pairs) {
  if (!optimotu_profile_enabled()) {
    return;
  }
  auto &t = tile_tls();
  t.active = true;
  t.index = index;
  t.pairs = pairs;
  t.start = std::chrono::steady_clock::now();
  t.have_pair_done = false;
  t.finish_s = 0;
}

void optimotu_profile_tile_mark_pair_done() {
  if (!optimotu_profile_enabled()) {
    return;
  }
  auto &t = tile_tls();
  if (!t.active || t.have_pair_done) {
    return;
  }
  t.pair_done = std::chrono::steady_clock::now();
  t.have_pair_done = true;
}

void optimotu_profile_tile_add_finish(double seconds) {
  if (!optimotu_profile_enabled()) {
    return;
  }
  tile_tls().finish_s += seconds;
}

void optimotu_profile_tile_leave() {
  if (!optimotu_profile_enabled()) {
    return;
  }
  auto &t = tile_tls();
  if (!t.active) {
    return;
  }
  const auto end = std::chrono::steady_clock::now();
  if (!t.have_pair_done) {
    t.pair_done = end;
    t.have_pair_done = true;
    optimotu_profile_mark_tile_work_done();
  }
  ProfileTile tile;
  tile.tile_index = t.index;
  tile.pairs = t.pairs;
  tile.pair_s = std::chrono::duration<double>(t.pair_done - t.start).count();
  tile.finish_s = t.finish_s;
  tile.lock_wait_s = optimotu_profile_take_lock_wait();
  optimotu_profile_take_which(tile.which_calls, tile.which_sum);
  optimotu_profile_take_applied(tile.applied);
  optimotu_profile_take_ksw2(tile.ksw2_aligned, tile.ksw2_skipped);
  t.active = false;
  optimotu_profile_add_tile(tile);
}

void optimotu_profile_report_merge(int n_seq, int threads) {
  if (!optimotu_profile_enabled()) {
    return;
  }
  auto &s = state();
  std::lock_guard<std::mutex> lock(s.mu);
  const double last_tile = s.last_tile_work_s.load();
  const double done = s.have_merge_complete ? s.merge_complete_s : 0;
  const double tail = (s.have_merge_complete && last_tile > 0)
                          ? std::max(0.0, done - last_tile)
                          : 0;
  double merge_cpu = 0;
  double merge_after_last_cpu = 0;
  double last_tile_merge_s = 0;
  for (const auto &m : s.merges) {
    merge_cpu += std::max(0.0, m.end_s - m.start_s);
    if (last_tile > 0 && m.end_s > last_tile) {
      merge_after_last_cpu += m.end_s - std::max(m.start_s, last_tile);
    }
    if (last_tile > 0 && m.start_s >= last_tile - 1e-9) {
      last_tile_merge_s = std::max(
          last_tile_merge_s, m.end_s - m.start_s);
    }
  }
  std::ostringstream line;
  line << "OPTIMOTU_MERGE_PROFILE"
       << " n_seq=" << n_seq
       << " threads=" << threads
       << " n_tile_marks=" << s.n_tile_marks.load()
       << " n_merges=" << s.merges.size()
       << " last_tile_work_s=" << last_tile
       << " merge_complete_s=" << done
       << " merge_tail_s=" << tail
       << " merge_cpu_s=" << merge_cpu
       << " merge_after_last_cpu_s=" << merge_after_last_cpu
       << " last_tile_merge_s=" << last_tile_merge_s
       << " parent_finalize_s=" << s.finalize_s;
  OPTIMOTU_CERR << line.str() << std::endl;
  if (const char *path = std::getenv("OPTIMOTU_MERGE_PROFILE_FILE")) {
    if (path[0] != '\0') {
      std::ofstream out(path, std::ios::app);
      if (out) {
        out << line.str() << std::endl;
      }
    }
  }
}

void optimotu_profile_report(
    int n_seq,
    int n_subsets,
    int threads,
    double score_wall_s,
    double sum_n2) {
  if (!optimotu_profile_enabled()) {
    return;
  }
  auto &s = state();
  std::lock_guard<std::mutex> lock(s.mu);
  double setup = 0;
  double pairs_s = 0;
  double finish = 0;
  double lock_wait = 0;
  double dist_s = 0;
  double relevant_s = 0;
  double apply_s = 0;
  std::uint64_t which_calls = 0;
  std::uint64_t which_sum = 0;
  std::uint64_t pairs_at_min = 0;
  std::uint64_t pairs_at_max = 0;
  double wall_min = 0;
  double wall_max = 0;
  bool have_wall = false;
  double setup_max = 0;
  double pair_max = 0;
  double finish_max = 0;
  std::uint64_t pairs = 0;
  std::uint64_t applied = 0;
  std::uint64_t samples = 0;
  double tile_cpu = 0;
  std::vector<double> tile_wall;
  tile_wall.reserve(s.tiles.size());
  for (const auto &t : s.tiles) {
    setup += t.setup_s;
    pairs_s += t.pair_s;
    finish += t.finish_s;
    lock_wait += t.lock_wait_s;
    setup_max = std::max(setup_max, t.setup_s);
    pair_max = std::max(pair_max, t.pair_s);
    finish_max = std::max(finish_max, t.finish_s);
    dist_s += t.dist_sample_s;
    relevant_s += t.relevant_sample_s;
    apply_s += t.apply_sample_s;
    pairs += t.pairs;
    applied += t.applied;
    samples += t.samples;
    which_calls += t.which_calls;
    which_sum += t.which_sum;
    const double wall = t.setup_s + t.pair_s + t.finish_s;
    tile_cpu += wall;
    if (!have_wall || wall < wall_min) {
      wall_min = wall;
      pairs_at_min = t.pairs;
    }
    if (!have_wall || wall > wall_max) {
      wall_max = wall;
      pairs_at_max = t.pairs;
    }
    have_wall = true;
    tile_wall.push_back(wall);
  }
  std::sort(tile_wall.begin(), tile_wall.end());
  const double tile_max = tile_wall.empty() ? 0 : tile_wall.back();
  const double tile_min = tile_wall.empty() ? 0 : tile_wall.front();
  const double tile_p50 = percentile_sorted(tile_wall, 0.5);
  const double sample_den = dist_s + relevant_s + apply_s;
  const double dist_frac = sample_den > 0 ? dist_s / sample_den : 0;
  const double relevant_frac = sample_den > 0 ? relevant_s / sample_den : 0;
  const double apply_frac = sample_den > 0 ? apply_s / sample_den : 0;
  const double mean_which = which_calls > 0
      ? static_cast<double>(which_sum) / static_cast<double>(which_calls)
      : 0;
  const double util =
      (threads > 0 && s.worker_s > 0)
          ? tile_cpu / (static_cast<double>(threads) * s.worker_s)
          : 0;
  const double tile_ratio = tile_min > 0 ? tile_max / tile_min : 0;
  OPTIMOTU_CERR
      << "OPTIMOTU_PROFILE"
      << " n_seq=" << n_seq
      << " n_subsets=" << n_subsets
      << " threads=" << threads
      << " tiles=" << s.tiles.size()
      << " worker_wall_s=" << s.worker_s
      << " finalize_wall_s=" << s.finalize_s
      << " score_wall_s=" << score_wall_s
      << " setup_cpu_s=" << setup
      << " pair_cpu_s=" << pairs_s
      << " finish_cpu_s=" << finish
      << " lock_wait_s=" << lock_wait
      << " setup_max_s=" << setup_max
      << " pair_max_s=" << pair_max
      << " finish_max_s=" << finish_max
      << " pairs=" << pairs
      << " applied=" << applied
      << " dist_frac=" << dist_frac
      << " relevant_frac=" << relevant_frac
      << " apply_frac=" << apply_frac
      << " mean_which=" << mean_which
      << " dist_samples=" << samples
      << " pairs_at_tile_min=" << pairs_at_min
      << " pairs_at_tile_max=" << pairs_at_max
      << " tile_min_s=" << tile_min
      << " tile_p50_s=" << tile_p50
      << " tile_max_s=" << tile_max
      << " tile_max_over_min=" << tile_ratio
      << " util=" << util
      << " sum_n2=" << sum_n2
      << std::endl;
  if (const char *path = std::getenv("OPTIMOTU_PROFILE_FILE")) {
    if (path[0] != '\0') {
      std::ofstream out(path, std::ios::app);
      if (out) {
        out << "OPTIMOTU_PROFILE"
            << " n_seq=" << n_seq
            << " n_subsets=" << n_subsets
            << " threads=" << threads
            << " tiles=" << s.tiles.size()
            << " worker_wall_s=" << s.worker_s
            << " finalize_wall_s=" << s.finalize_s
            << " pair_cpu_s=" << pairs_s
            << " finish_cpu_s=" << finish
            << " applied=" << applied
            << " mean_which=" << mean_which
            << " tile_min_s=" << tile_min
            << " tile_p50_s=" << tile_p50
            << " tile_max_s=" << tile_max
            << " tile_max_over_min=" << tile_ratio
            << " util=" << util
            << std::endl;
      }
    }
  }
  write_tile_tsv(s.tiles);
}
