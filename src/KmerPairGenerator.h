#ifndef OPTIMOTU_KMERPAIRGENERATOR_H_INCLUDED
#define OPTIMOTU_KMERPAIRGENERATOR_H_INCLUDED

#include <vector>
#include <unordered_map>
#include <memory>
#include <cstdint>

#include "PairGenerator.h"

// uses a kmer (8-mer) index over a set of sequences to return only the sequence
// pairs which are likely to be similar
class KmerPairGenerator : public PairGenerator {
protected:
  const std::shared_ptr<std::vector<std::vector<std::size_t>>> kmer_seq_index;
  const std::shared_ptr<std::vector<std::vector<std::uint16_t>>> seq_kmer_index;
  const double udist_threshold;
  std::size_t i;
  std::unordered_map<std::size_t, std::uint16_t> match_index; // matches for seq[i]
  std::unordered_map<std::size_t, std::uint16_t>::iterator j; // iterator through match_index

  // generate match_index() for seq[i]; i if there are no matches for seq[i]
  // then increment i until there are matches or end is reached.
  void update_match_index();

  // protected constructor for use in subset() which reuses the kmer index
  KmerPairGenerator(
    const std::size_t begin,
    const std::size_t end,
    const std::shared_ptr<std::vector<std::vector<std::size_t>>> kmer_seq_index,
    const std::shared_ptr<std::vector<std::vector<std::uint16_t>>> seq_kmer_index,
    const double udist_threshold
  );

public:
  // public constructor generates the kmer index
  KmerPairGenerator(
    const std::size_t begin,
    const std::size_t end,
    const std::vector<std::string> & seq,
    const double udist_threshold
  );

  std::unique_ptr<PairGenerator> subset(
      const std::size_t begin,
      const std::size_t end
  ) const override;

  bool operator()(std::pair<std::size_t, std::size_t> & pair) override;

};

#endif // OPTIMOTU_KMERPAIRGENERATOR_H_INCLUDED
