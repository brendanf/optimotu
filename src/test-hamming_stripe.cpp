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
