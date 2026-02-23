#include "optimotu.h"
#include "BipartitePairGenerator.h"
#include <cmath>

BipartitePairGenerator& BipartitePairGenerator::operator++() {
  if (++_j >= nj) {
    if (++_i >= ni + nj) {
      has_more = false;
    }
    _j = 0;
  }
  return *this;
}

std::size_t BipartitePairGenerator::forward_map(const std::size_t value) const {
  if (value < nj) {
    return value + begin_j;
  } else if (value < ni + nj) {
    return value - nj + begin_i;
  } else {
    OPTIMOTU_STOP(
      "Invalid value: " + std::to_string(value) +
      " for BipartitePairGenerator::forward_map() with ni = " + std::to_string(ni) +
      " and nj = " + std::to_string(nj)
    );
  }
}

std::size_t BipartitePairGenerator::reverse_map(const std::size_t value) const {
  if (value >= begin_i && value < end_i) {
    return value - begin_i + nj;
  } else if (value >= begin_j && value < end_j) {
    return value - begin_j;
  } else {
    OPTIMOTU_STOP(
      "Invalid value: " + std::to_string(value) +
      " for BipartitePairGenerator::reverse_map() with ni = " + std::to_string(ni) +
      " and nj = " + std::to_string(nj)
    );
  }
}

std::size_t BipartitePairGenerator::i0() const {
  return _i - nj + begin_i;
}
std::size_t BipartitePairGenerator::j0() const {
  return _j + begin_j;
}

#include <testthat.h>
#ifdef TESTTHAT_ENABLED
#include <set>

context("BipartitePairGenerator") {

  test_that("BipartitePairGenerator(0,0,0,0) has no pairs") {
    BipartitePairGenerator gen(0, 0, 0, 0);
    expect_false(gen);
  }

  test_that("BipartitePairGenerator 1x1 produces one pair") {
    BipartitePairGenerator gen(10, 11, 20, 21);
    expect_true(gen);
    expect_true(gen.i0() == 10);
    expect_true(gen.j0() == 20);
    ++gen;
    expect_false(gen);
  }

  test_that("BipartitePairGenerator(0,2,10,13) produces 6 pairs") {
    BipartitePairGenerator gen(0, 2, 10, 13);
    std::set<std::pair<std::size_t, std::size_t>> pairs;
    while (gen) {
      pairs.insert({gen.i0(), gen.j0()});
      ++gen;
    }
    expect_true(pairs.size() == 6u);
    for (std::size_t i = 0; i < 2; ++i) {
      for (std::size_t j = 10; j < 13; ++j) {
        expect_true(pairs.count({i, j}));
      }
    }
  }

  test_that("BipartitePairGenerator forward_map and reverse_map round-trip") {
    BipartitePairGenerator gen(5, 8, 1, 4);
    for (std::size_t v = 0; v <= gen.max_value(); ++v) {
      expect_true(gen.reverse_map(gen.forward_map(v)) == v);
    }
  }
}
#endif