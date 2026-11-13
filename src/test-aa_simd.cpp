#include <testthat.h>
#include "PackedAaSequenceSet.h"

#include <cmath>
#include <cstdint>
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

  test_that("evo column scan matches the scalar pair")
  {
    const char *letters = "ARNDCQEGHILKMFPSTWYV";
    EvoModel model;
    model.t_max = 5.0;
    model.t_grid = {0.0, 0.1, 0.25, 0.5, 1.0, 2.0, 5.0};
    model.lambda.assign(EvoModel::kState, -1.0);
    model.evec.assign(EvoModel::kState * EvoModel::kState, 0.0);
    model.evec_inv.assign(EvoModel::kState * EvoModel::kState, 0.0);
    for (int i = 0; i < EvoModel::kState; ++i)
    {
      model.lambda[static_cast<std::size_t>(i)] = -static_cast<double>(i + 1);
      model.evec[static_cast<std::size_t>(i) * EvoModel::kState +
                 static_cast<std::size_t>(i)] = 1.0;
      model.evec_inv[static_cast<std::size_t>(i) * EvoModel::kState +
                     static_cast<std::size_t>(i)] = 1.0;
      const int code = letters[i] - 'A' + 1;
      model.code_to_state[static_cast<std::size_t>(code)] = i;
      model.is_state[static_cast<std::size_t>(code)] = 1;
    }
    const int n_grid = model.n_grid();
    model.logp.assign(
        static_cast<std::size_t>(n_grid) * EvoModel::kCode * EvoModel::kCode,
        0.0f);
    for (int g = 0; g < n_grid; ++g)
    {
      const double t = model.t_grid[static_cast<std::size_t>(g)];
      for (int i = 0; i < EvoModel::kState; ++i)
      {
        const int ca = letters[i] - 'A' + 1;
        for (int j = 0; j < EvoModel::kState; ++j)
        {
          const int cb = letters[j] - 'A' + 1;
          const double ll = (i == j)
                                ? -t
                                : -((t - 0.25) * (t - 0.25) + 0.1);
          model.logp[(
                         (static_cast<std::size_t>(ca) * EvoModel::kCode +
                          static_cast<std::size_t>(cb)) *
                         static_cast<std::size_t>(n_grid)) +
                     static_cast<std::size_t>(g)] = static_cast<float>(ll);
        }
      }
    }
    std::vector<std::string> seq(20);
    uint32_t s = 7;
    for (int i = 0; i < 20; ++i)
    {
      seq[static_cast<std::size_t>(i)].assign(48, 'A');
      for (int k = 0; k < 48; ++k)
      {
        s = s * 1664525u + 1013904223u;
        seq[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] =
            letters[(s >> 27) % 20];
      }
    }
    seq[1] = std::string(8, '-') + seq[1].substr(8);
    PackedAaSequenceSet pss(seq);
    std::vector<int> idx(seq.size());
    for (int i = 0; i < static_cast<int>(seq.size()); ++i)
    {
      idx[static_cast<std::size_t>(i)] = i;
    }
    for (int q = 0; q < pss.num_seqs; ++q)
    {
      std::vector<uint8_t> pos;
      std::vector<int> targets;
      for (int j = 0; j < pss.num_seqs; ++j)
      {
        if (j == q)
        {
          continue;
        }
        targets.push_back(j);
      }
      pss.fill_pos_major_idx(targets.data(), static_cast<int>(targets.size()), pos);
      std::vector<double> dist(targets.size(), 0.0);
      std::vector<int> n_ok(targets.size(), 0);
      pss.score_query_evo(
          q,
          pos.data(),
          targets.data(),
          static_cast<int>(targets.size()),
          static_cast<int>(targets.size()),
          0,
          true,
          model,
          dist.data(),
          n_ok.data());
      for (int t = 0; t < static_cast<int>(targets.size()); ++t)
      {
        const double ref = pss.dist_evo(
            q, targets[static_cast<std::size_t>(t)], 0, true, model);
        expect_true(std::fabs(dist[static_cast<std::size_t>(t)] - ref) < 1e-6);
      }
    }
  }
}
