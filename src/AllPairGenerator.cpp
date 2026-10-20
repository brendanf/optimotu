// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "optimotu.h"
#include "AllPairGenerator.h"
#include "BipartitePairGenerator.h"
#include "MappedClusterAlgorithm.h"
#include "tile_pairs.h"
#include <algorithm>
#include <cmath>

AllPairGenerator::Builder::Builder(
    const std::size_t n,
    const std::size_t min_n_subgenerators) : DivisiblePairGenerator::Builder(n, min_n_subgenerators)
{
  // Square tiling: n_blocks diagonal tiles plus two pieces per
  // rectangle gives n_blocks^2 equal-work generators.
  n_tiles = std::max<std::size_t>(
      1,
      static_cast<std::size_t>(std::ceil(std::sqrt(static_cast<double>(
          std::max<std::size_t>(min_n_subgenerators, 1))))));
  n_subgenerators = n_tiles * n_tiles;
}

namespace
{

  void push_rectangle(
      std::vector<std::unique_ptr<PairGenerator>> &generators,
      std::size_t begin_i,
      std::size_t end_i,
      std::size_t begin_j,
      std::size_t end_j)
  {
    if (end_i <= begin_i || end_j <= begin_j)
    {
      return;
    }
    generators.push_back(std::make_unique<BipartitePairGenerator>(
        begin_i, end_i, begin_j, end_j));
  }

  void push_split_rectangle(
      std::vector<std::unique_ptr<PairGenerator>> &generators,
      std::size_t begin_i,
      std::size_t end_i,
      std::size_t begin_j,
      std::size_t end_j)
  {
    const std::size_t ni = end_i - begin_i;
    const std::size_t nj = end_j - begin_j;
    if (ni == 0 || nj == 0)
    {
      return;
    }
    if (ni >= 2 && ni >= nj)
    {
      const std::size_t mid = begin_i + ni / 2;
      push_rectangle(generators, begin_i, mid, begin_j, end_j);
      push_rectangle(generators, mid, end_i, begin_j, end_j);
      return;
    }
    if (nj >= 2)
    {
      const std::size_t mid = begin_j + nj / 2;
      push_rectangle(generators, begin_i, end_i, begin_j, mid);
      push_rectangle(generators, begin_i, end_i, mid, end_j);
      return;
    }
    push_rectangle(generators, begin_i, end_i, begin_j, end_j);
  }

} // namespace

AllPairGenerator &AllPairGenerator::operator++()
{
  if (++_j >= _i)
  {
    ++_i;
    _j = 0;
  }
  if (_i >= n)
    has_more = false;
  return *this;
}

std::size_t AllPairGenerator::i0() const { return _i + offset; }
std::size_t AllPairGenerator::j0() const { return _j + offset; }

std::size_t AllPairGenerator::forward_map(const std::size_t value) const
{
  return value + offset;
}
std::size_t AllPairGenerator::reverse_map(const std::size_t value) const
{
  return value - offset;
}

std::vector<std::unique_ptr<PairGenerator>> AllPairGenerator::Builder::build(int verbose) const
{
  OPTIMOTU_VERBOSE(
      1,
      << "Building " << n_subgenerators
      << " subgenerators for AllPairGenerator with offset " << offset
      << std::endl;);
  OPTIMOTU_VERBOSE(
      4,
      << "  Main diagonal subgenerators:" << std::endl;);
  std::vector<std::unique_ptr<PairGenerator>> generators;
  generators.reserve(n_subgenerators);

  // Create the subgenerators along the main diagonal
  double range = n;
  for (std::size_t tile_i = 0; tile_i < n_tiles; tile_i++)
  {
    std::size_t begin_i = tile_i * range / n_tiles;
    std::size_t end_i = (tile_i + 1) * range / n_tiles;
    if (end_i > begin_i)
    {
      OPTIMOTU_VERBOSE(
          4,
          << "    Subgenerator " << tile_i
          << ": " << begin_i << " to " << end_i
          << " (size = " << end_i - begin_i
          << ", offset = " << offset + begin_i
          << ")"
          << std::endl;);
      generators.push_back(
          std::make_unique<AllPairGenerator>(end_i - begin_i, offset + begin_i));
    }
  }
  OPTIMOTU_VERBOSE(
      4,
      << "  Subgenerators below the diagonal:" << std::endl;);
  // Create the subgenerators below the diagonal
  for (std::size_t tile_i = 1; tile_i < n_tiles; tile_i++)
  {
    std::size_t begin_i = tile_i * range / n_tiles;
    std::size_t end_i = (tile_i + 1) * range / n_tiles;
    if (end_i == begin_i)
      continue;
    for (std::size_t tile_j = 0; tile_j < tile_i; tile_j++)
    {
      std::size_t begin_j = tile_j * range / n_tiles;
      std::size_t end_j = (tile_j + 1) * range / n_tiles;
      if (end_j == begin_j)
        continue;
      OPTIMOTU_VERBOSE(
          4,
          << "    Subgenerator " << tile_i << "," << tile_j
          << ": " << begin_i << " (" << offset + begin_i
          << ") to " << end_i << " (" << offset + end_i
          << ") and " << begin_j << " (" << offset + begin_j
          << ") to " << end_j << " (" << offset + end_j
          << "), sizes = " << end_i - begin_i
          << ", " << end_j - begin_j
          << std::endl;);
      push_split_rectangle(
          generators,
          offset + begin_i,
          offset + end_i,
          offset + begin_j,
          offset + end_j);
    }
  }
  return generators;
}

#include <testthat.h>
#ifdef TESTTHAT_ENABLED
#include <set>

// Spec: pairs (i, j) with 1 <= i < n and 0 <= j < i (strict lower triangle).
context("AllPairGenerator") {

  test_that("AllPairGenerator(0) has no pairs") {
    AllPairGenerator gen(0);
    expect_false(gen);
  }

  test_that("AllPairGenerator(1) has no pairs") {
    AllPairGenerator gen(1);
    expect_false(gen);
  }

  test_that("AllPairGenerator(2) produces (1,0) then exhausts") {
    AllPairGenerator gen(2);
    expect_true(gen);
    expect_true(gen.i() == 1);
    expect_true(gen.j() == 0);
    expect_true(gen.i0() == 1);
    expect_true(gen.j0() == 0);
    ++gen;
    expect_false(gen);
  }

  test_that("AllPairGenerator(4) produces exactly 6 pairs") {
    AllPairGenerator gen(4);
    std::set<std::pair<std::size_t, std::size_t>> pairs;
    while (gen) {
      pairs.insert({gen.i0(), gen.j0()});
      ++gen;
    }
    expect_true(pairs.size() == 6u);
    expect_true(pairs.count({1, 0}));
    expect_true(pairs.count({2, 0}));
    expect_true(pairs.count({2, 1}));
    expect_true(pairs.count({3, 0}));
    expect_true(pairs.count({3, 1}));
    expect_true(pairs.count({3, 2}));
  }

  test_that("AllPairGenerator with offset") {
    const std::size_t offset = 10;
    AllPairGenerator gen(3, offset);
    expect_true(gen.i0() == offset + 1);
    expect_true(gen.j0() == offset);
    std::set<std::pair<std::size_t, std::size_t>> pairs;
    while (gen) {
      pairs.insert({gen.i0(), gen.j0()});
      ++gen;
    }
    expect_true(pairs.size() == 3u);
    expect_true(gen.forward_map(0) == offset);
    expect_true(gen.reverse_map(offset) == 0);
  }

  test_that("AllPairGenerator forward_map and reverse_map round-trip") {
    AllPairGenerator gen(5, 3);
    for (std::size_t v = 0; v < gen.max_value() + 1; ++v) {
      expect_true(gen.reverse_map(gen.forward_map(v)) == v);
    }
  }

  test_that("AllPairGenerator::Builder(10) subgenerators cover 45 pairs") {
    AllPairGenerator::Builder b(10, 4);
    b.set_offset(0);
    auto gens = b.build(0);
    std::set<std::pair<std::size_t, std::size_t>> all_pairs;
    for (auto& g : gens) {
      while (*g) {
        all_pairs.insert({g->i0(), g->j0()});
        ++(*g);
      }
    }
    expect_true(all_pairs.size() == 45u);
    for (std::size_t i = 1; i < 10; ++i) {
      for (std::size_t j = 0; j < i; ++j) {
        expect_true(all_pairs.count({i, j}));
      }
    }
  }
}

context("subset_tile_fwd_map")
{

  test_that("diagonal tile keeps sorted parent-locals in intersection")
  {
    std::unordered_map<j_t, j_t> global_to_subset{
        {0, 0}, {2, 1}, {4, 2}, {6, 3}, {8, 4}};
    AllPairGenerator tile(5, 0);
    auto intersection_fwd = subset_tile_fwd_map(global_to_subset, tile);
    expect_true(intersection_fwd.size() == 3u);
    expect_true(intersection_fwd[0] == 0u);
    expect_true(intersection_fwd[1] == 1u);
    expect_true(intersection_fwd[2] == 2u);
  }

  test_that("bipartite tile keeps sorted parent-locals in intersection")
  {
    std::unordered_map<j_t, j_t> global_to_subset{
        {0, 0}, {2, 1}, {4, 2}, {6, 3}, {8, 4}};
    BipartitePairGenerator tile(5, 10, 0, 5);
    auto intersection_fwd = subset_tile_fwd_map(global_to_subset, tile);
    expect_true(intersection_fwd.size() == 5u);
    expect_true(intersection_fwd[0] == 0u);
    expect_true(intersection_fwd[1] == 1u);
    expect_true(intersection_fwd[2] == 2u);
    expect_true(intersection_fwd[3] == 3u);
    expect_true(intersection_fwd[4] == 4u);
  }
}

context("DivisiblePairGenerator::Builder") {

  test_that("Builder respects min_n_subgenerators") {
    using B = AllPairGenerator::Builder;
    B b1(100, 1);
    expect_true(b1.build(0).size() >= 1u);
    B b5(100, 5);
    expect_true(b5.build(0).size() >= 5u);
  }
}

context("for_each_tile_pair") {

  test_that("visitor over Builder tiles covers the lower triangle") {
    for (std::size_t n : {2u, 5u, 10u}) {
      for (std::size_t min_gen : {1u, 4u, 8u}) {
        AllPairGenerator::Builder b(n, min_gen);
        auto gens = b.build(0);
        std::set<std::pair<std::size_t, std::size_t>> from_visitor;
        for (auto &g : gens) {
          for_each_tile_pair(*g, [&](
            std::size_t i0, std::size_t j0, std::size_t, std::size_t
          ) {
            from_visitor.insert({i0, j0});
          });
        }
        expect_true(from_visitor.size() == n * (n - 1) / 2);
        for (std::size_t i = 1; i < n; ++i) {
          for (std::size_t j = 0; j < i; ++j) {
            expect_true(from_visitor.count({i, j}));
          }
        }
      }
    }
  }

  test_that("for_each_assigned_tile walks a tile range without overlap")
  {
    AllPairGenerator::Builder b(10, 4);
    auto gens = b.build(0);
    std::vector<std::size_t> seen(gens.size(), 0);
    for (std::size_t begin = 0; begin < gens.size(); ++begin)
    {
      const std::size_t end = begin + 1;
      for_each_assigned_tile(
          gens, begin, end,
          [&](std::size_t pg_index, PairGenerator &)
          {
            ++seen[pg_index];
          });
    }
    for (std::size_t c : seen) {
      expect_true(c == 1u);
    }
  }

  test_that("Builder tiles have similar pair counts")
  {
    AllPairGenerator::Builder b(40, 9);
    auto gens = b.build(0);
    expect_true(gens.size() >= 9u);
    std::vector<std::size_t> counts;
    counts.reserve(gens.size());
    for (auto &g : gens)
    {
      std::size_t n_pairs = 0;
      for_each_tile_pair(*g, [&](
                                 std::size_t, std::size_t, std::size_t, std::size_t)
                         { ++n_pairs; });
      if (n_pairs > 0)
      {
        counts.push_back(n_pairs);
      }
    }
    expect_true(!counts.empty());
    std::sort(counts.begin(), counts.end());
    const std::size_t median = counts[counts.size() / 2];
    expect_true(median > 0);
    expect_true(counts.back() <= median * 3);
    expect_true(counts.front() * 3 >= median);
  }
}
#endif