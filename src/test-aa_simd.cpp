#include <testthat.h>
#include "PackedAaSequenceSet.h"

#include <cmath>
#include <string>
#include <vector>

namespace {

void expect_near(double a, double b) {
  expect_true(std::fabs(a - b) < 1e-12);
}

int identity_site_ref(uint8_t a, uint8_t b, bool ignore_gaps) {
  if (a == 0 && b == 0) {
    return 0;
  }
  if (ignore_gaps && (a == 0 || b == 0)) {
    return 0;
  }
  return 1;
}

}  // namespace

context("aa_simd") {
  test_that("column identity matches the scalar loop") {
    std::vector<std::string> seq = {
      "ACDEFGHIKLMNPQRSTVWYACDE",
      "ACDEFGHIKLMNPQRSTVWYACDE",
      "ACDEFGHIKLMNPQRSTVWY----",
      "-CDEFGHIKLMNPQRSTVWYACDE",
      "AAAAAAAAAAAAAAAAAAAAAAAA"
    };
    PackedAaSequenceSet pss(seq);
    std::vector<uint8_t> pos;
    pss.fill_pos_major(1, pss.num_seqs, pos);
    const int n_t = pss.num_seqs - 1;
    std::vector<double> dist(static_cast<std::size_t>(n_t), 1.0);
    std::vector<int> n_ok(static_cast<std::size_t>(n_t), 0);
    std::vector<int> idx(static_cast<std::size_t>(n_t));
    for (int t = 0; t < n_t; ++t) {
      idx[static_cast<std::size_t>(t)] = t + 1;
    }
    pss.score_query_vs_pos_major(
      0, pos.data(), idx.data(), n_t, n_t, 0, true, nullptr, 1.0, false,
      dist.data(), n_ok.data()
    );
    for (int t = 0; t < n_t; ++t) {
      expect_near(dist[static_cast<std::size_t>(t)], pss.dist(0, t + 1, 0, true, nullptr));
    }
  }

  test_that("target counts not a multiple of 16 still match") {
    std::vector<std::string> seq(19, std::string(17, 'A'));
    seq[3][8] = 'C';
    seq[18][0] = '-';
    PackedAaSequenceSet pss(seq);
    std::vector<uint8_t> pos;
    pss.fill_pos_major(1, 19, pos);
    std::vector<double> dist(18, 1.0);
    std::vector<int> n_ok(18, 0);
    std::vector<int> idx(18);
    for (int t = 0; t < 18; ++t) {
      idx[static_cast<std::size_t>(t)] = t + 1;
    }
    pss.score_query_vs_pos_major(
      0, pos.data(), idx.data(), 18, 18, 0, true, nullptr, 1.0, false,
      dist.data(), n_ok.data()
    );
    for (int t = 0; t < 18; ++t) {
      expect_near(dist[static_cast<std::size_t>(t)], pss.dist(0, t + 1, 0, true, nullptr));
    }
    (void)identity_site_ref;
  }

  test_that("ignore_gaps false with unequal ranges matches pairwise") {
    std::vector<std::string> seq = {
      "--DEFGHIKLMNPQRSTVWYACDE",
      "ACDEFGHIKLMNPQRSTVWY----",
      "ACDEFGHIKLMNPQRSTVWYACDE",
      "-CDEFGHIK-----STVWYACDE-",
      "AAAAAAAAAAAAAAAAAAAAAAAA"
    };
    PackedAaSequenceSet pss(seq);
    std::vector<uint8_t> pos;
    pss.fill_pos_major(1, pss.num_seqs, pos);
    const int n_t = pss.num_seqs - 1;
    std::vector<double> dist(static_cast<std::size_t>(n_t), 1.0);
    std::vector<int> n_ok(static_cast<std::size_t>(n_t), 0);
    std::vector<int> idx(static_cast<std::size_t>(n_t));
    for (int t = 0; t < n_t; ++t) {
      idx[static_cast<std::size_t>(t)] = t + 1;
    }
    pss.score_query_vs_pos_major(
      0, pos.data(), idx.data(), n_t, n_t, 0, false, nullptr, 1.0, false,
      dist.data(), n_ok.data()
    );
    for (int t = 0; t < n_t; ++t) {
      expect_near(
        dist[static_cast<std::size_t>(t)],
        pss.dist(0, t + 1, 0, false, nullptr)
      );
    }

    const int n = PackedAaSequenceSet::kLut;
    std::vector<float> lut(static_cast<std::size_t>(n * n), 1.0f);
    for (int i = 0; i < n; ++i) {
      lut[static_cast<std::size_t>(i) * n + static_cast<std::size_t>(i)] = 0.0f;
    }
    dist.assign(static_cast<std::size_t>(n_t), 1.0);
    n_ok.assign(static_cast<std::size_t>(n_t), 0);
    pss.score_query_vs_pos_major(
      0, pos.data(), idx.data(), n_t, n_t, 0, false, lut.data(), 1.0, false,
      dist.data(), n_ok.data()
    );
    for (int t = 0; t < n_t; ++t) {
      expect_near(
        dist[static_cast<std::size_t>(t)],
        pss.dist(0, t + 1, 0, false, lut.data())
      );
    }
  }
}
