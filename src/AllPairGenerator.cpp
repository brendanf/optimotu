// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "optimotu.h"
#include "AllPairGenerator.h"
#include <cmath>

AllPairGenerator& AllPairGenerator::operator++() {
  if (++_j >= _i) {
    ++_i;
    _j = 0;
  }
  if (_i >= n) has_more = false;
  return *this;
}

std::size_t AllPairGenerator::i0() const { return _i + offset; }
std::size_t AllPairGenerator::j0() const { return _j + offset; }

std::size_t AllPairGenerator::forward_map(const std::size_t value) const {
  return value + offset;
}
std::size_t AllPairGenerator::reverse_map(const std::size_t value) const {
  return value - offset;
}

std::vector<std::unique_ptr<PairGenerator>> AllPairGenerator::Builder::build(int verbose) const {
  OPTIMOTU_VERBOSE(
    1,
    << "Building " << n_subgenerators
    << " subgenerators for AllPairGenerator with offset " << offset
    << std::endl;
  );
  OPTIMOTU_VERBOSE(
    4,
    << "  Main diagonal subgenerators:" << std::endl;
  );
  std::vector<std::unique_ptr<PairGenerator>> generators;
  generators.reserve(n_subgenerators);

  // Create the subgenerators along the main diagonal
  double range = n;
  for (std::size_t tile_i = 0; tile_i < n_tiles; tile_i++) {
    std::size_t begin_i = tile_i * range / n_tiles;
    std::size_t end_i = (tile_i + 1) * range / n_tiles;
    if (end_i > begin_i) {
      OPTIMOTU_VERBOSE(
        4,
        << "    Subgenerator " << tile_i
        << ": " << begin_i << " to " << end_i
        << " (size = " << end_i - begin_i
        << ", offset = " << offset + begin_i
        << ")"
        << std::endl;
      );
      generators.push_back(
        std::make_unique<AllPairGenerator>(end_i - begin_i, offset + begin_i)
      );
    }
  }
  OPTIMOTU_VERBOSE(
    4,
    << "  Subgenerators below the diagonal:" << std::endl;
  );
  // Create the subgenerators below the diagonal
  for (std::size_t tile_i = 1; tile_i < n_tiles; tile_i++) {
    std::size_t begin_i = tile_i * range / n_tiles;
    std::size_t end_i = (tile_i + 1) * range / n_tiles;
    if (end_i == begin_i) continue;
    for (std::size_t tile_j = 0; tile_j < tile_i; tile_j++) {
      std::size_t begin_j = tile_j * range / n_tiles;
      std::size_t end_j = (tile_j + 1) * range / n_tiles;
      if (end_j == begin_j) continue;
      OPTIMOTU_VERBOSE(
        4,
        << "    Subgenerator " << tile_i << "," << tile_j
        << ": " << begin_i << " (" << offset + begin_i
        << ") to " << end_i << " (" << offset + end_i
        << ") and " << begin_j << " (" << offset + begin_j
        << ") to " << end_j << " (" << offset + end_j
        << "), sizes = " << end_i - begin_i
        << ", " << end_j - begin_j
        << std::endl;
      );
      generators.push_back(
        std::make_unique<BipartitePairGenerator>(
          offset + begin_i,
          offset + end_i,
          offset + begin_j,
          offset + end_j
        )
      );
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

context("DivisiblePairGenerator::Builder") {

  test_that("Builder respects min_n_subgenerators") {
    using B = AllPairGenerator::Builder;
    B b1(100, 1);
    expect_true(b1.build(0).size() >= 1u);
    B b5(100, 5);
    expect_true(b5.build(0).size() >= 5u);
  }
}
#endif