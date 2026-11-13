// SPDX-FileCopyrightText: 2026 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_PROFILE_STATS_H_INCLUDED
#define OPTIMOTU_PROFILE_STATS_H_INCLUDED

#include <cstddef>
#include <cstdint>

// OPTIMOTU_PROFILE=1 prints one phase line per clustering batch.
// OPTIMOTU_PROFILE_TILES=<path> appends one TSV row per tile.
// The hot path checks this once and skips the extra clocks otherwise.

bool optimotu_profile_enabled();

struct ProfileTile {
  std::size_t tile_index = 0;
  double setup_s = 0;
  double pair_s = 0;
  double finish_s = 0;
  double lock_wait_s = 0;
  double dist_sample_s = 0;
  double relevant_sample_s = 0;
  double apply_sample_s = 0;
  std::uint64_t pairs = 0;
  std::uint64_t applied = 0;
  std::uint64_t samples = 0;
  std::uint64_t which_calls = 0;
  std::uint64_t which_sum = 0;
  std::uint64_t ksw2_aligned = 0;
  std::uint64_t ksw2_skipped = 0;
};

void optimotu_profile_reset();
void optimotu_profile_add_tile(const ProfileTile &tile);
// Thread-local. Called on the worker that holds the merge-cache lock.
void optimotu_profile_add_lock_wait(double seconds);
double optimotu_profile_take_lock_wait();
void optimotu_profile_add_which(std::uint32_t n_subsets);
void optimotu_profile_take_which(
    std::uint64_t &calls,
    std::uint64_t &sum);
void optimotu_profile_add_applied(std::uint64_t n);
void optimotu_profile_take_applied(std::uint64_t &n);
void optimotu_profile_add_ksw2(std::uint64_t aligned, std::uint64_t skipped);
void optimotu_profile_take_ksw2(
    std::uint64_t &aligned,
    std::uint64_t &skipped);
void optimotu_profile_set_worker(double seconds);
void optimotu_profile_set_finalize(double seconds);
// Wall clock for the merge tail: last tile's local work vs merge done.
void optimotu_profile_mark_origin();
double optimotu_profile_now_s();
void optimotu_profile_mark_tile_work_done();
void optimotu_profile_add_merge_span(double start_s, double end_s);
void optimotu_profile_mark_merge_complete();
void optimotu_profile_tile_enter(std::size_t index, std::uint64_t pairs);
void optimotu_profile_tile_mark_pair_done();
void optimotu_profile_tile_add_finish(double seconds);
void optimotu_profile_tile_leave();
void optimotu_profile_report_merge(int n_seq, int threads);
void optimotu_profile_report(
    int n_seq,
    int n_subsets,
    int threads,
    double score_wall_s,
    double sum_n2);

#endif
