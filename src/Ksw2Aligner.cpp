// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#include "Ksw2Aligner.h"

#include <algorithm>
#include <cstring>
#include <sstream>

namespace {

char cigar_op_char(uint32_t op) {
  switch (op) {
  case KSW_CIGAR_EQ:
    return '=';
  case KSW_CIGAR_X:
    return 'X';
  case KSW_CIGAR_MATCH:
    return 'M';
  case KSW_CIGAR_INS:
    return 'I';
  case KSW_CIGAR_DEL:
    return 'D';
  case KSW_CIGAR_N_SKIP:
    return 'N';
  default:
    return '?';
  }
}

} // namespace

Ksw2Aligner::Ksw2Aligner(
  int match,
  int mismatch,
  int gap_open,
  int gap_extend,
  int gap_open2,
  int gap_extend2
) : match(match),
    mismatch(mismatch),
    gap_open(gap_open),
    gap_extend(gap_extend),
    gap_open2(gap_open2),
    gap_extend2(gap_extend2) {
  std::memset(&ez, 0, sizeof(ez));
  rebuild_matrix();
}

Ksw2Aligner::~Ksw2Aligner() {
  free_cigar();
}

void Ksw2Aligner::free_cigar() {
  if (ez.cigar) {
    kfree(0, ez.cigar);
    ez.cigar = nullptr;
  }
  ez.n_cigar = 0;
  ez.m_cigar = 0;
}

void Ksw2Aligner::rebuild_matrix() {
  // KSW2 needs a>0 and b<0. R match==0 (default) maps to +1 so edit-like
  // defaults behave like Edlib / WFA2 under identity distance.
  sc_match = static_cast<int8_t>(match > 0 ? match : 1);
  sc_mismatch = static_cast<int8_t>(
    mismatch > 0 ? -mismatch : (mismatch < 0 ? mismatch : -1)
  );
  // Mirror WFA2 GapAffine vs GapAffine2Pieces selection: piece-2 disabled
  // when both open2 and extend2 are zero; identical pieces collapse to
  // single affine.
  dual_affine =
    !((gap_open2 == 0 && gap_extend2 == 0) ||
      (gap_open == gap_open2 && gap_extend == gap_extend2));

  // 5x5 DNA matrix; last residue is wildcard (N / other).
  const int m = 5;
  for (int i = 0; i < m - 1; ++i) {
    for (int j = 0; j < m - 1; ++j) {
      mat[i * m + j] = (i == j) ? sc_match : sc_mismatch;
    }
    mat[i * m + (m - 1)] = 0;
    mat[(m - 1) * m + i] = 0;
  }
  mat[(m - 1) * m + (m - 1)] = 0;
}

void Ksw2Aligner::encode(
  const char *seq,
  int len,
  std::vector<uint8_t> &out
) const {
  out.resize(static_cast<size_t>(len));
  for (int i = 0; i < len; ++i) {
    switch (seq[i]) {
    case 'A':
    case 'a':
      out[i] = 0;
      break;
    case 'C':
    case 'c':
      out[i] = 1;
      break;
    case 'G':
    case 'g':
      out[i] = 2;
      break;
    case 'T':
    case 't':
    case 'U':
    case 'u':
      out[i] = 3;
      break;
    default:
      out[i] = 4;
      break;
    }
  }
}

bool Ksw2Aligner::align(
  const char *a,
  int a_len,
  const char *b,
  int b_len,
  AlignmentSpan span
) {
  free_cigar();
  ksw_reset_extz(&ez);
  last_ok = false;
  last_score = 0;

  if (a_len <= 0 || b_len <= 0) {
    return false;
  }

  encode(a, a_len, query_enc);
  encode(b, b_len, target_enc);
  last_query = query_enc;
  last_target = target_enc;

  const int8_t q = static_cast<int8_t>(gap_open);
  const int8_t e = static_cast<int8_t>(gap_extend);
  const int8_t q2 = static_cast<int8_t>(gap_open2);
  const int8_t e2 = static_cast<int8_t>(gap_extend2);
  const int w = band_width;
  const int zdrop = -1;
  const int end_bonus = 0;
  int flag = 0;
  if (span == AlignmentSpan::EXTEND) {
    flag |= KSW_EZ_EXTZ_ONLY;
  }

#ifdef __SSE2__
  if (dual_affine) {
    ksw_extd2_sse(
      0, a_len, query_enc.data(), b_len, target_enc.data(),
      5, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, &ez
    );
  } else {
    ksw_extz2_sse(
      0, a_len, query_enc.data(), b_len, target_enc.data(),
      5, mat, q, e, w, zdrop, end_bonus, flag, &ez
    );
  }
#else
  if (dual_affine) {
    ksw_extd(
      0, a_len, query_enc.data(), b_len, target_enc.data(),
      5, mat, q, e, q2, e2, w, zdrop, flag, &ez
    );
  } else {
    ksw_extz(
      0, a_len, query_enc.data(), b_len, target_enc.data(),
      5, mat, q, e, w, zdrop, flag, &ez
    );
  }
#endif

  if (ez.zdropped) {
    return false;
  }

  // Convert MATCH runs to EQ/X using the encoded sequences. Do not use
  // ksw_cigar2eqx(): it does not update *ci1 after krealloc.
  if (ez.n_cigar > 0 && ez.cigar) {
    int nc1 = 0;
    int mc1 = 0;
    uint32_t *ci1 = nullptr;
    int x = 0;
    int y = 0;
    for (int k = 0; k < ez.n_cigar; ++k) {
      const int op = ez.cigar[k] & 0xf;
      const int len = static_cast<int>(ez.cigar[k] >> 4);
      if (op == KSW_CIGAR_MATCH) {
        for (int i = 0; i < len; ++i) {
          const uint32_t eqx =
            (target_enc[static_cast<size_t>(x + i)] ==
             query_enc[static_cast<size_t>(y + i)])
              ? KSW_CIGAR_EQ
              : KSW_CIGAR_X;
          ci1 = ksw_push_cigar(0, &nc1, &mc1, ci1, eqx, 1);
        }
        x += len;
        y += len;
      } else {
        ci1 = ksw_push_cigar(0, &nc1, &mc1, ci1, op, len);
        if (op == KSW_CIGAR_DEL || op == KSW_CIGAR_N_SKIP) {
          x += len;
        } else if (op == KSW_CIGAR_INS) {
          y += len;
        } else if (op == KSW_CIGAR_EQ || op == KSW_CIGAR_X) {
          x += len;
          y += len;
        }
      }
    }
    kfree(0, ez.cigar);
    ez.cigar = ci1;
    ez.n_cigar = nc1;
    ez.m_cigar = mc1;
  }

  if (span == AlignmentSpan::EXTEND) {
    // Prefer end-reaching extension score when available.
    last_score = ez.reach_end ? ez.mqe : static_cast<int>(ez.max);
  } else {
    last_score = ez.score;
  }
  if (last_score == KSW_NEG_INF) {
    return false;
  }
  if (ez.n_cigar <= 0 || ez.cigar == nullptr) {
    return false;
  }
  last_ok = true;
  return true;
}

std::string Ksw2Aligner::getCIGAR() const {
  if (!last_ok || ez.n_cigar <= 0 || ez.cigar == nullptr) {
    return "";
  }
  std::ostringstream ss;
  for (int i = 0; i < ez.n_cigar; ++i) {
    const uint32_t c = ez.cigar[i];
    const int len = static_cast<int>(c >> 4);
    const char op = cigar_op_char(c & 0xf);
    if (len > 1) {
      ss << len;
    }
    ss << op;
  }
  return ss.str();
}
