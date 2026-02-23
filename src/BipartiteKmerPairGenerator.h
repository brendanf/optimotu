#ifndef OPTIMOTU_BIPARTITEKMERPAIRGENERATOR_H_INCLUDED
#define OPTIMOTU_BIPARTITEKMERPAIRGENERATOR_H_INCLUDED

#include "PairGenerator.h"
#include <cstdint>
#include <memory>
#include <vector>
#include <unordered_map>

using kmer_seq_index_t = std::vector<std::vector<std::size_t>>;
using seq_kmer_index_t = std::vector<std::vector<std::uint16_t>>;

class BipartiteKmerPairGenerator : public PairGenerator {
  // offsets for the indices
  const std::size_t begin_i, end_i, begin_j, end_j, ni, nj;
  // threshold for considering a kmer to be similar
  const double udist_threshold;
  // sequence-to-kmer index: for each sequence in [begin_i, end_i), which kmers
  // are found in it?
  const std::shared_ptr<seq_kmer_index_t> seq_kmer_index_i;
  // sequence-to-kmer index: for each sequence in [begin_j, end_j), which kmers
  // are found in it?
  const std::shared_ptr<seq_kmer_index_t> seq_kmer_index_j;
  // kmer-to-sequence index: for each kmer, which of the sequences in [begin_j, end_j)
  // are they found in?
  const std::shared_ptr<kmer_seq_index_t> kmer_seq_index_j;
  // matches for seq[i]: how many times does each sequence match seq[i]?
  std::unordered_map<std::size_t, std::uint16_t> match_index;
  // iterator through match_index
  std::unordered_map<std::size_t, std::uint16_t>::iterator match_index_j;
  // generate match_index for seq[i]; if there are no matches for seq[i]
  // then increment i until there are matches or end is reached.
  void update_match_index();

public:
  // This constructor accepts the kmer-to-sequence and sequence-to-kmer indices
  // They should be pre-computed to cover only the sequences [begin_j, end_j)
  // and [begin_i, end_i) respectively, and re-indexed to the mapped range.
  // (This is not done by the constructor, to allow re-use of the indices for
  // multiple generators.)
  BipartiteKmerPairGenerator(
    const std::size_t begin_i,
    const std::size_t end_i,
    const std::size_t begin_j,
    const std::size_t end_j,
    const double udist_threshold,
    const std::shared_ptr<seq_kmer_index_t> seq_kmer_index_i,
    const std::shared_ptr<seq_kmer_index_t> seq_kmer_index_j,
    const std::shared_ptr<kmer_seq_index_t> kmer_seq_index_j
  ) :
  PairGenerator(end_i - begin_i + end_j - begin_j),
  begin_i(begin_i), end_i(end_i),
  begin_j(begin_j), end_j(end_j),
  ni(end_i - begin_i), nj(end_j - begin_j),
  udist_threshold(udist_threshold),
  seq_kmer_index_i(seq_kmer_index_i),
  seq_kmer_index_j(seq_kmer_index_j),
  kmer_seq_index_j(kmer_seq_index_j) {
    _i = nj;  // internal row index: [nj, ni+nj) are i-side
    has_more = (ni > 0 && nj > 0);
    update_match_index();
    if (match_index.empty()) has_more = false;
  };

  BipartiteKmerPairGenerator& operator++() override;

  std::size_t forward_map(const std::size_t value) const override;

  std::size_t reverse_map(const std::size_t value) const override;

  std::size_t i0() const override;
  std::size_t j0() const override;

  std::size_t max_value() const override;

};

#endif // OPTIMOTU_BIPARTITEKMERPAIRGENERATOR_H_INCLUDED