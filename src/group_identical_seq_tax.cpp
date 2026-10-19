// SPDX-CopyrightText: (c) 2025, Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifdef OPTIMOTU_R

#include <Rcpp.h>
#include <cstdint>
#include <cstring>
#include <unordered_map>
#include <vector>

namespace {

uint64_t fnv1a(const char *p, std::size_t n) {
  uint64_t h = 14695981039346656037ull;
  for (std::size_t i = 0; i < n; ++i) {
    h ^= static_cast<unsigned char>(p[i]);
    h *= 1099511628211ull;
  }
  return h;
}

uint64_t mix_hash(uint64_t h, uint64_t x) {
  h ^= x + 0x9e3779b97f4a7c15ull + (h << 6) + (h >> 2);
  return h;
}

bool charsxp_equal(SEXP a, SEXP b) {
  if (a == b) return true;
  if (a == NA_STRING || b == NA_STRING) return false;
  const R_xlen_t na = LENGTH(a);
  const R_xlen_t nb = LENGTH(b);
  if (na != nb) return false;
  return std::memcmp(CHAR(a), CHAR(b), static_cast<std::size_t>(na)) == 0;
}

} // namespace

//' Group byte-identical sequences that share all rank labels
//'
//' @param seq (`character`) reference sequences
//' @param ranks (`data.frame`) rank columns aligned with `seq`
//' @return a list with `representative` (1-based original indices of
//'   first-seen group members) and `map` (1-based index into
//'   `representative` for every original row)
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List group_identical_seq_tax(
    Rcpp::CharacterVector seq,
    Rcpp::DataFrame ranks
) {
  const R_xlen_t n = seq.size();
  if (ranks.nrow() != n) {
    Rcpp::stop(
      "rank rows (%d) must match sequence length (%d)",
      ranks.nrow(),
      n
    );
  }
  const int n_rank = ranks.size();
  std::vector<Rcpp::CharacterVector> rank_cols;
  rank_cols.reserve(static_cast<std::size_t>(n_rank));
  for (int r = 0; r < n_rank; ++r) {
    rank_cols.emplace_back(ranks[r]);
    if (rank_cols.back().size() != n) {
      Rcpp::stop("rank column length does not match sequences");
    }
  }

  auto row_hash = [&](R_xlen_t i) -> uint64_t {
    SEXP s = seq[i];
    uint64_t h = 0;
    if (s != NA_STRING) {
      h = fnv1a(CHAR(s), static_cast<std::size_t>(LENGTH(s)));
    }
    for (int r = 0; r < n_rank; ++r) {
      SEXP rnk = rank_cols[r][i];
      h = mix_hash(h, reinterpret_cast<uintptr_t>(rnk));
    }
    return h;
  };

  auto row_equal = [&](R_xlen_t a, R_xlen_t b) -> bool {
    if (!charsxp_equal(seq[a], seq[b])) return false;
    for (int r = 0; r < n_rank; ++r) {
      if (!charsxp_equal(rank_cols[r][a], rank_cols[r][b])) return false;
    }
    return true;
  };

  std::unordered_map<uint64_t, std::vector<int>> buckets;
  buckets.reserve(static_cast<std::size_t>(n));
  Rcpp::IntegerVector map(n);
  std::vector<int> representative;
  representative.reserve(static_cast<std::size_t>(n));

  for (R_xlen_t i = 0; i < n; ++i) {
    if (seq[i] == NA_STRING) {
      Rcpp::stop("sequences cannot contain NA");
    }
    const uint64_t h = row_hash(i);
    auto &bucket = buckets[h];
    int found = -1;
    for (int rep : bucket) {
      if (row_equal(i, rep)) {
        found = rep;
        break;
      }
    }
    if (found < 0) {
      bucket.push_back(static_cast<int>(i));
      representative.push_back(static_cast<int>(i + 1));
      map[i] = static_cast<int>(representative.size());
    } else {
      map[i] = map[found];
    }
  }

  Rcpp::IntegerVector rep_out(representative.begin(), representative.end());
  return Rcpp::List::create(
      Rcpp::Named("representative") = rep_out,
      Rcpp::Named("map") = map
  );
}

//' Remap subset membership onto sequence representatives
//'
//' @param seq_index (`list` of `integer`) 1-based original row indices
//' @param true_partition (`list` of `integer`) true labels parallel to
//'   `seq_index`
//' @param map (`integer`) 1-based representative index per original row
//' @param seq_names (`character`) names of the compacted sequences
//' @return a list with remapped `seq_index`, `true_partition`,
//'   `seq_id`, `multiplicity`, and `n_cluster`
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::List remap_testset_representatives(
    Rcpp::List seq_index,
    Rcpp::List true_partition,
    Rcpp::IntegerVector map,
    Rcpp::CharacterVector seq_names
) {
  const int n_sub = seq_index.size();
  if (true_partition.size() != n_sub) {
    Rcpp::stop("true_partition and seq_index must have the same length");
  }
  const int n_rep = seq_names.size();
  Rcpp::List out_index(n_sub);
  Rcpp::List out_part(n_sub);
  Rcpp::List out_id(n_sub);
  Rcpp::List out_mult(n_sub);
  Rcpp::IntegerVector n_cluster(n_sub);

  std::vector<int> counts;
  std::vector<int> first_part;
  std::vector<int> seen_at;
  counts.assign(static_cast<std::size_t>(n_rep), 0);
  first_part.assign(static_cast<std::size_t>(n_rep), NA_INTEGER);
  seen_at.assign(static_cast<std::size_t>(n_rep), -1);
  int stamp = 0;

  for (int s = 0; s < n_sub; ++s) {
    Rcpp::IntegerVector old_idx = seq_index[s];
    Rcpp::IntegerVector old_part = true_partition[s];
    const int n_old = old_idx.size();
    if (old_part.size() != n_old) {
      Rcpp::stop(
        "true_partition[[%d]] length does not match seq_index",
        s + 1
      );
    }
    ++stamp;
    std::vector<int> order;
    order.reserve(static_cast<std::size_t>(n_old));
    for (int j = 0; j < n_old; ++j) {
      const int orig = old_idx[j];
      if (orig == NA_INTEGER || orig < 1 || orig > map.size()) {
        Rcpp::stop("seq_index contains an out-of-range index");
      }
      const int m = map[orig - 1];
      if (m == NA_INTEGER || m < 1 || m > n_rep) {
        Rcpp::stop("representative map is out of range");
      }
      const std::size_t m0 = static_cast<std::size_t>(m - 1);
      if (seen_at[m0] != stamp) {
        seen_at[m0] = stamp;
        counts[m0] = 0;
        first_part[m0] = old_part[j];
        order.push_back(m);
      }
      counts[m0] += 1;
    }
    const int n_keep = static_cast<int>(order.size());
    Rcpp::IntegerVector new_idx(n_keep);
    Rcpp::IntegerVector new_part(n_keep);
    Rcpp::IntegerVector new_mult(n_keep);
    Rcpp::CharacterVector new_id(n_keep);
    for (int k = 0; k < n_keep; ++k) {
      const int m = order[static_cast<std::size_t>(k)];
      const std::size_t m0 = static_cast<std::size_t>(m - 1);
      new_idx[k] = m;
      new_part[k] = first_part[m0];
      new_mult[k] = counts[m0];
      SET_STRING_ELT(new_id, k, static_cast<SEXP>(seq_names[m - 1]));
    }
    out_index[s] = new_idx;
    out_part[s] = new_part;
    out_id[s] = new_id;
    out_mult[s] = new_mult;
    n_cluster[s] = n_keep;
  }

  return Rcpp::List::create(
      Rcpp::Named("seq_index") = out_index,
      Rcpp::Named("true_partition") = out_part,
      Rcpp::Named("seq_id") = out_id,
      Rcpp::Named("multiplicity") = out_mult,
      Rcpp::Named("n_cluster") = n_cluster
  );
}

#endif
