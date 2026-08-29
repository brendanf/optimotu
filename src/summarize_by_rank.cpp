// SPDX-CopyrightText: (c) 2025, Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifdef OPTIMOTU_R

#include <Rcpp.h>
#include <cstring>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <cstdint>

// R CHARSXPs are typically interned; identical strings share a pointer, so
// we can store IDs without copying bytes into std::string. Map keys still
// compare by character content so output row order matches the historical
// std::map<std::string, ...> behaviour.

struct CharsexpLess
{
  bool operator()(SEXP a, SEXP b) const
  {
    return std::strcmp(CHAR(a), CHAR(b)) < 0;
  }
};

// a single identification of a sequence as a taxon
struct TaxonID {
  SEXP seq_id;
  SEXP taxon;
  int seq_index; // 1-based row in input data
  TaxonID(SEXP seq_id, SEXP taxon, int seq_index)
      : seq_id(seq_id), taxon(taxon), seq_index(seq_index) {}
};

typedef std::vector<TaxonID> TaxonIDVec;

struct Subtaxa
{
  // Unique taxa; true_partition IDs follow lexicographic order (historical
  // std::set<std::string> behaviour).
  std::set<SEXP, CharsexpLess> taxa;
  TaxonIDVec seq_map;
  Subtaxa() : taxa(), seq_map() {}
};

typedef std::map<SEXP, Subtaxa, CharsexpLess> SupertaxonMap;

//' Summarize taxonomic ranks by superordinate rank
//' @param data ('data.frame') the taxonomy to summarize; should contain a
//' column named `seq_id` and columns for each value of `ranks`
//' @param ranks ('character') the ranks to summarize
//' @param id_col ('character') the name of the column in `data` containing
//'  sequence IDs; defaults to `seq_id`. Note that the *output* will always
//'  name the column `seq_id`.
//' @return a data frame with columns:
//'
//'  - `supertaxon` (`character`) the superordinate taxon
//'  - `superrank` (`character`) the rank of the superordinate taxon
//'  - `rank` (`character`) the rank being summarized
//'  - `n_taxa` (`integer`) the number of unique taxa at the rank
//'  - `n_seq` (`integer`) the number of sequences
//'  - `seq_id` (`list` of `character`) a list of sequence IDs
//'  - `seq_index` (`list` of `integer`) 1-based row indices in `data` for
//'     each element in `seq_id`
//'  - `true_partition` (`list` of `integer`) integer mapping to taxa for each
//'     element in `seq_id`
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::RObject summarize_by_rank(
    Rcpp::DataFrame data,
    Rcpp::CharacterVector ranks,
    std::string id_col = "seq_id")
{

  auto id_col_cstr = id_col.c_str();
  if (!data.containsElementNamed(id_col_cstr))
  {
    Rcpp::stop("data must contain a column named '%s'", id_col_cstr);
  }
  Rcpp::CharacterVector seq_id = data[id_col];
  if (Rcpp::any(Rcpp::is_na(seq_id)))
  {
    Rcpp::stop("data$%s cannot contain NA values", id_col_cstr);
  }

  std::vector<Rcpp::CharacterVector> taxonomy(ranks.size());
  for (int i = 0; i < ranks.size(); i++)
  {
    if (!data.containsElementNamed(ranks[i]))
    {
      Rcpp::String rank_i = ranks[i];
      Rcpp::stop("data must contain a column named '%s'", rank_i.get_cstring());
    }
    Rcpp::String rank_i = ranks[i];
    taxonomy[i] = data[rank_i];
  }

  // First vector index is the superordinate rank,
  // second vector index is the subordinate rank,
  // third map index is the superordinate taxon.
  std::vector<std::vector<SupertaxonMap>> taxon_map;

  for (int super_r = 0; super_r < ranks.size() - 1; super_r++)
  {
    taxon_map.emplace_back();
    for (int sub_r = super_r + 1; sub_r < ranks.size(); sub_r++)
    {
      taxon_map[super_r].emplace_back();
      auto &sub_map = taxon_map[super_r][sub_r - super_r - 1];
      for (int i = 0; i < data.nrow(); i++) {
        SEXP super_taxon = taxonomy[super_r][i];
        SEXP sub_taxon = taxonomy[sub_r][i];
        if (super_taxon == NA_STRING) continue;
        if (sub_taxon == NA_STRING) continue;

        auto map_entry = sub_map.try_emplace(super_taxon);
        Subtaxa &st = map_entry.first->second;
        st.taxa.insert(sub_taxon);
        st.seq_map.emplace_back(seq_id[i], sub_taxon, i + 1);
      }
    }
  }

  int n_out = 0;
  for (const auto & super_map : taxon_map) {
    for (const auto & sub_map : super_map) {
      n_out += static_cast<int>(sub_map.size());
    }
  }

  Rcpp::CharacterVector out_supertaxon(n_out);
  Rcpp::CharacterVector out_superrank(n_out);
  Rcpp::CharacterVector out_rank(n_out);
  Rcpp::IntegerVector n_taxa(n_out);
  Rcpp::IntegerVector n_seq(n_out);
  Rcpp::List seq_ids(n_out);
  Rcpp::List seq_indices(n_out);
  Rcpp::List true_taxa(n_out);

  int i = 0;
  for (int super_r = 0; super_r < ranks.size() - 1; super_r++) {
    for (int sub_r = super_r + 1; sub_r < ranks.size(); sub_r++) {
      const auto & sub_map = taxon_map[super_r][sub_r - super_r - 1];
      for (const auto &entry : sub_map)
      {
        SEXP super_taxon = entry.first;
        const Subtaxa &seq_taxa = entry.second;
        // Assign CHARSXP without allocating a new string.
        SET_STRING_ELT(out_supertaxon, i, super_taxon);
        SET_STRING_ELT(out_superrank, i, ranks[super_r]);
        SET_STRING_ELT(out_rank, i, ranks[sub_r]);
        n_taxa[i] = static_cast<int>(seq_taxa.taxa.size());
        n_seq[i] = static_cast<int>(seq_taxa.seq_map.size());
        Rcpp::CharacterVector seq_id_vec(seq_taxa.seq_map.size());
        Rcpp::IntegerVector seq_index_vec(seq_taxa.seq_map.size());
        Rcpp::IntegerVector true_taxa_vec(seq_taxa.seq_map.size());
        std::unordered_map<SEXP, int> taxon_id_map;
        taxon_id_map.reserve(seq_taxa.taxa.size());
        int tax_j = 0;
        for (SEXP taxon : seq_taxa.taxa)
        {
          taxon_id_map.emplace(taxon, tax_j);
          tax_j++;
        }
        for (std::size_t j = 0; j < seq_taxa.seq_map.size(); j++) {
          SET_STRING_ELT(seq_id_vec, j, seq_taxa.seq_map[j].seq_id);
          seq_index_vec[j] = seq_taxa.seq_map[j].seq_index;
          true_taxa_vec[j] = taxon_id_map.at(seq_taxa.seq_map[j].taxon);
        }
        seq_ids[i] = seq_id_vec;
        seq_indices[i] = seq_index_vec;
        true_taxa[i] = true_taxa_vec;
        i++;
      }
    }
  }

  Rcpp::List output = Rcpp::List::create(
      Rcpp::Named("supertaxon") = out_supertaxon,
      Rcpp::Named("superrank") = out_superrank,
      Rcpp::Named("rank") = out_rank,
      Rcpp::Named("n_taxa") = n_taxa,
      Rcpp::Named("n_seq") = n_seq,
      Rcpp::Named("seq_id") = seq_ids,
      Rcpp::Named("seq_index") = seq_indices,
      Rcpp::Named("true_partition") = true_taxa);

  output.attr("class") = Rcpp::CharacterVector::create("tbl_df", "tbl", "data.frame");
  output.attr("row.names") = Rcpp::seq_len(n_out);

  return output;
}

#endif
