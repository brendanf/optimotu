#ifndef OPTIMOTU_KMERPAIRGENERATOR_H_INCLUDED
#define OPTIMOTU_KMERPAIRGENERATOR_H_INCLUDED

#include <vector>
#include <unordered_map>
#include <memory>
#include <cstdint>
#include <string>

#include "PairGenerator.h"
#include "BipartiteKmerPairGenerator.h"

// uses a kmer (8-mer) index over a set of sequences to return only the sequence
// pairs which are likely to be similar
class KmerPairGenerator : public DivisiblePairGenerator {
protected:

  // kmer-to-sequence index: for each kmer, which sequences is it found in?
  const std::shared_ptr<kmer_seq_index_t> kmer_seq_index;
  // sequence-to-kmer index: for each sequence, which kmers are found in it?
  const std::shared_ptr<seq_kmer_index_t> seq_kmer_index;
  // matches for seq[i]: how many times does each sequence match seq[i]?
  std::unordered_map<std::size_t, std::uint16_t> match_index;
  // iterator through match_index
  std::unordered_map<std::size_t, std::uint16_t>::iterator match_index_j;
  // threshold for considering a match to be similar
  const double udist_threshold;
  // offset for the indices
  const std::size_t offset = 0;

  // generate match_index() for seq[i]; i if there are no matches for seq[i]
  // then increment i until there are matches or end is reached.
  void update_match_index();

public:
  // basic constructor generates the kmer index
  KmerPairGenerator(
    const std::vector<std::string> & seq,
    const double udist_threshold,
    const std::size_t offset = 0
  );

  // constructor for use in Builder::build() which reuses the kmer index
  KmerPairGenerator(
    const std::shared_ptr<kmer_seq_index_t> kmer_seq_index,
    const std::shared_ptr<seq_kmer_index_t> seq_kmer_index,
    const double udist_threshold,
    const std::size_t offset = 0
  );

  KmerPairGenerator& operator++() override;
  std::size_t forward_map(const std::size_t value) const override;
  std::size_t reverse_map(const std::size_t value) const override;
  std::size_t i0() const override;
  std::size_t j0() const override;

  class Builder : public DivisiblePairGenerator::Builder {
  protected:
    const std::vector<std::string> & seq;
    double udist_threshold;
    std::size_t offset;
  public:
    Builder(const std::vector<std::string> & seq, const std::size_t min_n_subgenerators = 1) :
      DivisiblePairGenerator::Builder(seq.size(), min_n_subgenerators), seq(seq) {};

    void set_udist_threshold(const double udist_threshold) { this->udist_threshold = udist_threshold; }
    void set_offset(const std::size_t offset) { this->offset = offset; }

    std::vector<std::unique_ptr<PairGenerator>> build(int verbose = 0) const override;
  };
};

#endif // OPTIMOTU_KMERPAIRGENERATOR_H_INCLUDED
