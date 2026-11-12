#include <testthat.h>
#include "PackedSequenceSet.h"
#include "hamming_stripe.h"
extern "C" {
#include "defs.h"
}

#include <cmath>
#include <string>
#include <vector>

namespace {

std::vector<std::string> dna_set(std::initializer_list<const char *> xs) {
  std::vector<std::string> out;
  out.reserve(xs.size());
  for (const char *s : xs) {
    out.emplace_back(s);
  }
  return out;
}

void expect_near(double a, double b) {
  expect_true(std::fabs(a - b) < 1e-12);
}

}  // namespace

context("hamming_stripe") {
  test_that("stripe matches pdistB for terminal gaps and N") {
    auto seq = dna_set({
      "--ACGTACGTNNNNACGTACGTACGTACGTACGT",
      "GTACGTACGTACGTACGTACGTACGTACGTAC--",
      "ACGTACGTACGTACGTACGTACGTACGTACGTAC",
      "NNNNACGTACGTACGTACGTACGTACGTACGTNN"
    });
    PackedSequenceSet pss(seq);
    std::vector<uint64_t> seq_major;
    std::vector<uint64_t> mask_major;
    pss.fill_word_major(1, pss.num_seqs, seq_major, mask_major);
    const int n_t = pss.num_seqs - 1;
    std::vector<double> dist(static_cast<std::size_t>(n_t));
    std::vector<int> n_ok(static_cast<std::size_t>(n_t));
    std::vector<int> n_match(static_cast<std::size_t>(n_t));
    pss.score_query_vs_word_major(
      0, seq_major.data(), mask_major.data(), n_t, n_t, 0,
      dist.data(), n_ok.data(), n_match.data()
    );
    for (int j = 1; j < pss.num_seqs; ++j) {
      int s = pss.start[0] >= pss.start[j] ? pss.start[0] : pss.start[j];
      int e = pss.end[0] <= pss.end[j] ? pss.end[0] : pss.end[j];
      int num_ok = 0;
      int num_matches = 0;
      double d = 1.0;
      if (s < e) {
        pdistB(
          pss.packed_seq[0].data(), pss.mask[0].data(),
          pss.packed_seq[j].data(), pss.mask[j].data(),
          s, e, 0, &num_ok, &num_matches, &d
        );
      }
      expect_true(n_ok[j - 1] == num_ok);
      expect_true(n_match[j - 1] == num_matches);
      expect_near(dist[j - 1], d);
    }
  }

  test_that("min_overlap failure yields distance 1") {
    auto seq = dna_set({"ACGTACGT", "ACGTNNNN"});
    PackedSequenceSet pss(seq);
    std::vector<uint64_t> seq_major;
    std::vector<uint64_t> mask_major;
    pss.fill_word_major(1, 2, seq_major, mask_major);
    double d = 0.0;
    int n_ok = 0;
    pss.score_query_vs_word_major(
      0, seq_major.data(), mask_major.data(), 1, 1, 8,
      &d, &n_ok, nullptr
    );
    expect_near(d, 1.0);
  }
}

context("hamming_stripe_tails") {
  test_that("length 21 and 5 targets match pdistB") {
    std::vector<std::string> seq(5, std::string(21, 'A'));
    seq[1][3] = 'C';
    seq[2][10] = 'G';
    seq[3][20] = 'T';
    seq[4][0] = 'N';
    PackedSequenceSet pss(seq);
    std::vector<uint64_t> seq_major;
    std::vector<uint64_t> mask_major;
    pss.fill_word_major(1, 5, seq_major, mask_major);
    std::vector<double> dist(4);
    std::vector<int> n_ok(4);
    std::vector<int> n_match(4);
    pss.score_query_vs_word_major(
      0, seq_major.data(), mask_major.data(), 4, 4, 0,
      dist.data(), n_ok.data(), n_match.data()
    );
    for (int j = 1; j < 5; ++j) {
      int s = pss.start[0] >= pss.start[j] ? pss.start[0] : pss.start[j];
      int e = pss.end[0] <= pss.end[j] ? pss.end[0] : pss.end[j];
      int num_ok = 0;
      int num_matches = 0;
      double d = 1.0;
      pdistB(
        pss.packed_seq[0].data(), pss.mask[0].data(),
        pss.packed_seq[j].data(), pss.mask[j].data(),
        s, e, 0, &num_ok, &num_matches, &d
      );
      expect_true(n_ok[j - 1] == num_ok);
      expect_true(n_match[j - 1] == num_matches);
      expect_near(dist[j - 1], d);
    }
  }

  test_that("swap-remove scores the reduced set") {
    std::vector<std::string> seq(4, std::string(32, 'A'));
    seq[2][4] = 'C';
    PackedSequenceSet pss(seq);
    std::vector<int> idx = {1, 2, 3};
    std::vector<uint64_t> seq_major;
    std::vector<uint64_t> mask_major;
    pss.fill_word_major_idx(idx.data(), 3, seq_major, mask_major);
    const int last = 2;
    for (int w = 0; w < pss.ulen; ++w) {
      std::swap(
        seq_major[static_cast<std::size_t>(w) * 3],
        seq_major[static_cast<std::size_t>(w) * 3 + last]
      );
    }
    for (int w = 0; w < pss.mulen; ++w) {
      std::swap(
        mask_major[static_cast<std::size_t>(w) * 3],
        mask_major[static_cast<std::size_t>(w) * 3 + last]
      );
    }
    std::swap(idx[0], idx[last]);
    const int live = 2;
    std::vector<double> dist(static_cast<std::size_t>(live));
    std::vector<int> n_ok(static_cast<std::size_t>(live));
    pss.score_query_vs_word_major(
      0, seq_major.data(), mask_major.data(), live, 3, 0,
      dist.data(), n_ok.data(), nullptr
    );
    for (int t = 0; t < live; ++t) {
      expect_near(dist[t], pss.dist(0, idx[t], 0, true));
    }
  }
}

context("hamming_stripe_internal") {
  auto check_pdistb2 = [](
    const PackedSequenceSet &pss,
    int min_overlap
  ) {
    std::vector<uint64_t> seq_major;
    std::vector<uint64_t> mask_major;
    pss.fill_word_major(1, pss.num_seqs, seq_major, mask_major);
    const int n_t = pss.num_seqs - 1;
    std::vector<double> dist(static_cast<std::size_t>(n_t));
    std::vector<int> n_ok(static_cast<std::size_t>(n_t));
    std::vector<int> n_match(static_cast<std::size_t>(n_t));
    pss.score_query_vs_word_major(
      0, seq_major.data(), mask_major.data(), n_t, n_t, min_overlap,
      dist.data(), n_ok.data(), n_match.data(),
      pss.start.data() + 1, pss.end.data() + 1
    );
    for (int j = 1; j < pss.num_seqs; ++j) {
      expect_near(dist[j - 1], pss.dist(0, j, min_overlap, false));
      int s = pss.start[0] >= pss.start[j] ? pss.start[0] : pss.start[j];
      int e = pss.end[0] <= pss.end[j] ? pss.end[0] : pss.end[j];
      int num_ok = 0;
      int num_matches = 0;
      double d = 1.0;
      if (s < e) {
        pdistB2(
          pss.packed_seq[0].data(), pss.mask[0].data(),
          pss.packed_seq[j].data(), pss.mask[j].data(),
          s, e, min_overlap, &num_ok, &num_matches, &d
        );
      }
      expect_true(n_ok[j - 1] == num_ok);
      expect_true(n_match[j - 1] == num_matches);
      expect_near(dist[j - 1], d);
    }
  };

  test_that("terminal gaps and internal N match pdistB2") {
    auto seq = dna_set({
      "--ACGTACGTNNNNACGTACGTACGTACGTACGT",
      "GTACGTACGTACGTACGTACGTACGTACGTAC--",
      "ACGTACGTACGTACGTACGTACGTACGTACGTAC",
      "NNNNACGTACGTACGTACGTACGTACGTACGTNN"
    });
    PackedSequenceSet pss(seq);
    check_pdistb2(pss, 0);
  }

  test_that("overlap shorter than 64 bases matches pdistB2") {
    auto seq = dna_set({
      "ACGTACGTACGTACGT",
      "--GTACGTACGTACGT",
      "ACGTACGTACGT----",
      "NNNNACGTACGTACGT"
    });
    PackedSequenceSet pss(seq);
    check_pdistb2(pss, 0);
  }

  test_that("window crossing a 64-base boundary matches pdistB2") {
    std::string a(80, 'A');
    std::string b(80, 'A');
    std::string c(80, 'A');
    for (int i = 0; i < 10; ++i) {
      a[static_cast<std::size_t>(i)] = '-';
    }
    for (int i = 70; i < 80; ++i) {
      b[static_cast<std::size_t>(i)] = '-';
    }
    c[50] = 'N';
    c[65] = 'C';
    PackedSequenceSet pss(dna_set({a.c_str(), b.c_str(), c.c_str()}));
    check_pdistb2(pss, 0);
  }

  test_that("end on a 64-base boundary matches pdistB2") {
    std::vector<std::string> seq(5, std::string(64, 'A'));
    seq[1][0] = '-';
    seq[2][63] = '-';
    seq[3][32] = 'N';
    seq[4][16] = 'C';
    PackedSequenceSet pss(seq);
    check_pdistb2(pss, 0);
  }

  test_that("multiples of 64 bases keep the final mask word")
  {
    auto check = [](int n, int mismatch_at)
    {
      std::vector<std::string> seq(
          2, std::string(static_cast<std::size_t>(n), 'A'));
      if (mismatch_at >= 0)
      {
        seq[1][static_cast<std::size_t>(mismatch_at)] = 'C';
      }
      PackedSequenceSet pss(seq);
      const int s = pss.start[0] >= pss.start[1] ? pss.start[0]
                                                 : pss.start[1];
      const int e = pss.end[0] <= pss.end[1] ? pss.end[0] : pss.end[1];
      int num_ok = 0;
      int num_matches = 0;
      double d = 1.0;
      pdistB2(
          pss.packed_seq[0].data(), pss.mask[0].data(),
          pss.packed_seq[1].data(), pss.mask[1].data(),
          s, e, 0, &num_ok, &num_matches, &d);
      const int expect_match = mismatch_at >= 0 ? n - 1 : n;
      expect_true(s == 0);
      expect_true(e == n);
      expect_true(num_ok == n);
      expect_true(num_matches == expect_match);
      expect_near(
          d, 1.0 - static_cast<double>(expect_match) / static_cast<double>(n));

      std::vector<uint64_t> seq_major;
      std::vector<uint64_t> mask_major;
      pss.fill_word_major(1, 2, seq_major, mask_major);
      double stripe_d = 0.0;
      int stripe_ok = 0;
      int stripe_match = 0;
      pss.score_query_vs_word_major(
          0, seq_major.data(), mask_major.data(), 1, 1, 0,
          &stripe_d, &stripe_ok, &stripe_match,
          pss.start.data() + 1, pss.end.data() + 1);
      expect_true(stripe_ok == num_ok);
      expect_true(stripe_match == num_matches);
      expect_near(stripe_d, d);
    };
    check(64, -1);
    check(64, 63);
    check(128, -1);
    check(128, 127);
    check(192, 100);
  }

  test_that("empty overlap yields distance 1") {
    auto seq = dna_set({
      "ACGTACGT----------------",
      "----------------ACGTACGT"
    });
    PackedSequenceSet pss(seq);
    check_pdistb2(pss, 0);
  }

  test_that("min_overlap failure yields distance 1") {
    auto seq = dna_set({"ACGTACGT", "ACGTNNNN"});
    PackedSequenceSet pss(seq);
    std::vector<uint64_t> seq_major;
    std::vector<uint64_t> mask_major;
    pss.fill_word_major(1, 2, seq_major, mask_major);
    double d = 0.0;
    int n_ok = 0;
    pss.score_query_vs_word_major(
      0, seq_major.data(), mask_major.data(), 1, 1, 8,
      &d, &n_ok, nullptr, pss.start.data() + 1, pss.end.data() + 1
    );
    expect_near(d, 1.0);
    expect_near(d, pss.dist(0, 1, 8, false));
  }
}
