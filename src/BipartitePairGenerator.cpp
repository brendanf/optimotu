#include "optimotu.h"
#include "BipartitePairGenerator.h"
#include <cmath>

BipartitePairGenerator::BipartitePairGenerator(
    const std::size_t begin_i,
    const std::size_t end_i,
    const std::size_t begin_j,
    const std::size_t end_j) : PairGenerator(end_i - begin_i + end_j - begin_j),
                               begin_i(begin_i),
                               end_i(end_i),
                               begin_j(begin_j),
                               end_j(end_j),
                               ni(end_i - begin_i),
                               nj(end_j - begin_j)
{
  // Internal _i runs from nj to ni+nj-1 (rows); _j from 0 to nj-1
  _i = nj;
  has_more = (ni > 0 && nj > 0);
  if (has_more && end_j > begin_i)
  {
    OPTIMOTU_STOP(
        "BipartitePairGenerator requires end_j <= begin_i so that j0 < i0 "
        "(begin_i=" +
        std::to_string(begin_i) + ", end_i=" + std::to_string(end_i) +
        ", begin_j=" + std::to_string(begin_j) + ", end_j=" +
        std::to_string(end_j) + ")");
  }
}

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
    BipartitePairGenerator gen(20, 21, 10, 11);
    expect_true(gen);
    expect_true(gen.i0() == 20);
    expect_true(gen.j0() == 10);
    ++gen;
    expect_false(gen);
  }

  test_that("BipartitePairGenerator(10,13,0,2) produces 6 pairs") {
    BipartitePairGenerator gen(10, 13, 0, 2);
    std::set<std::pair<std::size_t, std::size_t>> pairs;
    while (gen) {
      pairs.insert({gen.i0(), gen.j0()});
      ++gen;
    }
    expect_true(pairs.size() == 6u);
    for (std::size_t i = 10; i < 13; ++i) {
      for (std::size_t j = 0; j < 2; ++j) {
        expect_true(pairs.count({i, j}));
      }
    }
  }

  test_that("inverted BipartitePairGenerator stops") {
    expect_error(BipartitePairGenerator(0, 2, 10, 13));
  }

  test_that("BipartitePairGenerator forward_map and reverse_map round-trip") {
    BipartitePairGenerator gen(5, 8, 1, 4);
    for (std::size_t v = 0; v <= gen.max_value(); ++v) {
      expect_true(gen.reverse_map(gen.forward_map(v)) == v);
    }
  }
}
#endif