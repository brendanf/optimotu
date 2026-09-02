#include <testthat.h>
#include "wfa_identity_bound.h"
#include <cmath>

context("wfa_identity_bound") {
  test_that("edit penalties recover ceil(maxd1)+1") {
    const double l1 = 317.0;
    const double l2 = 317.0;
    const double tau = 0.4;
    const double e_max = tau * (l1 + l2) / (2.0 - tau);
    WfaIdentityBound b = wfa_identity_bound(
        l1, l2, tau, 0, 1, 0, 1, 0, 1
    );
    expect_true(b.max_alignment_steps == (int)std::ceil(e_max) + 1);
    const double sigma = 1.0 - tau;
    const double denom = 1.0 + sigma;
    expect_true(
        b.max_k == (int)std::ceil((l2 - l1 * sigma) / denom)
    );
    expect_true(
        b.min_k == -(int)std::ceil((l1 - l2 * sigma) / denom)
    );
  }

  test_that("single affine uses mismatch and gap vertices") {
    const double l1 = 317.0;
    const double l2 = 317.0;
    const double tau = 0.4;
    const int x = 6;
    const int o = 4;
    const int e = 2;
    const int unit = o + e;
    const double delta = 0.0;
    const double x_max = tau * l2 - delta;
    const double e_max = tau * (l1 + l2) / (2.0 - tau);
    const int s1 = x * (int)std::ceil(x_max) + unit * (int)std::ceil(delta);
    const int s2 = unit * (int)std::ceil(e_max);
    const int smax = s1 > s2 ? s1 : s2;
    WfaIdentityBound b = wfa_identity_bound(
        l1, l2, tau, 0, x, o, e, o, e
    );
    expect_true(b.max_alignment_steps == smax + 1);
    expect_true(s2 > s1);
  }

  test_that("length difference is a mandatory insertion cost") {
    const double l1 = 100.0;
    const double l2 = 120.0;
    const double tau = 0.4;
    WfaIdentityBound edit = wfa_identity_bound(
        l1, l2, tau, 0, 1, 0, 1, 0, 1
    );
    WfaIdentityBound aff = wfa_identity_bound(
        l1, l2, tau, 0, 6, 4, 2, 4, 2
    );
    expect_true(aff.max_alignment_steps > edit.max_alignment_steps);
    expect_true(aff.max_k == edit.max_k);
    expect_true(aff.min_k == edit.min_k);
  }

  test_that("zero threshold only allows exact matches plus Delta gaps") {
    WfaIdentityBound b = wfa_identity_bound(
        10.0, 10.0, 0.0, 0, 6, 4, 2, 4, 2
    );
    expect_true(b.max_alignment_steps == 1);
  }
}
