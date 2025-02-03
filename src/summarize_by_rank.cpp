// SPDX-CopyrightText: (c) 2025, Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifdef OPTIMOTU_R

#include <Rcpp.h>
#include <map>
#include <set>
#include <vector>
#include <cstdint>

// a single identification of a sequence as a taxon
struct TaxonID {
  std::string seq_id;
  std::string taxon;
  TaxonID(std::string seq_id, std::string taxon) : seq_id(seq_id), taxon(taxon) {}
};

// a list of sequence IDs and their taxonomic identification
typedef std::vector<TaxonID> TaxonIDVec;

// a list of unique taxa (presumed to be subtaxa of a given supertaxon) and
// sequence IDs which map to them
struct Subtaxa {
  std::set<std::string> taxa;
  TaxonIDVec seq_map;
  Subtaxa() : taxa(), seq_map() {}
};

// a map linking the superordinate taxon to the subtaxa
typedef std::map<std::string, Subtaxa> SupertaxonMap;

//' Summarize taxonomic ranks by superordinate rank
//' @param data ('data.frame') the taxonomy to summarize; should contain a
//' column named `seq_id` and columns for each value of `ranks`
//' @param ranks ('character') the ranks to summarize
//' @return a data frame with columns:
//'
//'  - `supertaxon` (`character`) the superordinate taxon
//'  - `superrank` (`character`) the rank of the superordinate taxon
//'  - `rank` (`character`) the rank being summarized
//'  - `n_taxa` (`integer`) the number of unique taxa at the rank
//'  - `n_seq` (`integer`) the number of sequences
//'  - `seq_id` (`list` of `character`) a list of sequence IDs
//'  - `true_parition` (`list` of `integer`) integer mapping to taxa for each
//'     element in `seq_id`
//' @keywords internal
//' @export
// [[Rcpp::export]]
Rcpp::RObject summarize_by_rank(
    Rcpp::DataFrame data,
    Rcpp::CharacterVector ranks
) {

  // access the columns by name
  if (!data.containsElementNamed("seq_id")) {
    Rcpp::stop("data must contain a column named 'seq_id'");
  }
  Rcpp::CharacterVector seq_id = data["seq_id"];
  if (Rcpp::any(Rcpp::is_na(seq_id))) {
    Rcpp::stop("data$seq_id cannot contain NA values");
  }
  std::vector<std::string> seq_id_str = Rcpp::as<std::vector<std::string>>(seq_id);

  std::vector<Rcpp::CharacterVector> taxonomy(ranks.size());
  for (int i = 0; i < ranks.size(); i++) {
    if (!data.containsElementNamed(ranks[i])) {
      Rcpp::String rank_i = ranks[i];
      Rcpp::stop("data must contain a column named '%s'", rank_i.get_cstring());
    }
    Rcpp::String rank_i = ranks[i];
    taxonomy[i] = data[rank_i];
  }

  // First vector index is the superordinate rank,
  // second vector index is the subordinate rank,
  // third map index is the superordinate taxon
  // fourth pair contains:
  //    set of unique taxa at subordinate rank,
  //    vector of pairs of sequence ID and taxon (at subordinate rank)
  std::vector<std::vector<SupertaxonMap>> taxon_map;

  // fill the taxonomy map
  for (int super_r = 0; super_r < ranks.size() - 1; super_r++) {
    taxon_map.emplace_back();
    for (int sub_r = super_r + 1; sub_r < ranks.size(); sub_r++) {
      taxon_map[super_r].emplace_back();
      for (int i = 0; i < data.nrow(); i++) {
        Rcpp::String super_taxon = taxonomy[super_r][i];
        Rcpp::String sub_taxon = taxonomy[sub_r][i];
        if (super_taxon == NA_STRING) continue;
        if (sub_taxon == NA_STRING) continue;

        std::string super_taxon_str(super_taxon.get_cstring());
        std::string sub_taxon_str(sub_taxon.get_cstring());

        auto map_entry = taxon_map[super_r][sub_r - super_r - 1].try_emplace(super_taxon_str);
        map_entry.first->second.taxa.insert(sub_taxon_str);
        map_entry.first->second.seq_map.emplace_back(seq_id_str[i], sub_taxon_str);
      }
    }
  }

  // count the size of the output data frame
  int n_out = 0;
  for (const auto & super_map : taxon_map) {
    for (const auto & sub_map : super_map) {
      n_out += sub_map.size();
    }
  }

  // allocate the output data frame
  Rcpp::CharacterVector supertaxon(n_out);
  Rcpp::CharacterVector superrank(n_out);
  Rcpp::CharacterVector rank(n_out);
  Rcpp::IntegerVector n_taxa(n_out);
  Rcpp::IntegerVector n_seq(n_out);
  Rcpp::List seq_ids(n_out);
  Rcpp::List true_taxa(n_out);

  // Fill the output
  int i = 0;
  for (int super_r = 0; super_r < ranks.size() - 1; super_r++) {
    for (int sub_r = super_r + 1; sub_r < ranks.size(); sub_r++) {
      const auto & sub_map = taxon_map[super_r][sub_r - super_r - 1];
      for (const auto & [super_taxon, seq_taxa] : sub_map) {
        supertaxon[i] = super_taxon;
        superrank[i] = ranks[super_r];
        rank[i] = ranks[sub_r];
        n_taxa[i] = seq_taxa.taxa.size();
        n_seq[i] = seq_taxa.seq_map.size();
        Rcpp::CharacterVector seq_id_vec(seq_taxa.seq_map.size());
        Rcpp::IntegerVector true_taxa_vec(seq_taxa.seq_map.size());
        std::map<std::string, int> taxon_id_map;
        int j = 0;
        for (const auto & taxon : seq_taxa.taxa) {
          taxon_id_map[taxon] = j;
          j++;
        }
        for (std::size_t j = 0; j < seq_taxa.seq_map.size(); j++) {
          seq_id_vec[j] = seq_taxa.seq_map[j].seq_id;
          true_taxa_vec[j] = taxon_id_map[seq_taxa.seq_map[j].taxon];
        }
        seq_ids[i] = seq_id_vec;
        true_taxa[i] = true_taxa_vec;
        i++;
      }
    }
  }

  // generate a List, because Rcpp::DataFrame does not supposrt list columns
  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("supertaxon") = supertaxon,
    Rcpp::Named("superrank") = superrank,
    Rcpp::Named("rank") = rank,
    Rcpp::Named("n_taxa") = n_taxa,
    Rcpp::Named("n_seq") = n_seq,
    Rcpp::Named("seq_id") = seq_ids,
    Rcpp::Named("true_partition") = true_taxa
  );

  // make it a valid tibble
  output.attr("class") = Rcpp::CharacterVector::create("tbl_df", "tbl", "data.frame");
  output.attr("row.names") = Rcpp::seq_len(n_out);

  return output;
}

#endif
