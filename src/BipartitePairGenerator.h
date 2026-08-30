#ifndef OPTIMOTU_BIPARTITEPAIRGENERATOR_H_INCLUDED
#define OPTIMOTU_BIPARTITEPAIRGENERATOR_H_INCLUDED

#include "PairGenerator.h"

// Generates all pairs {i, j} where nj <= i < ni + nj and 0 <= j < nj.
// The forward mapping function maps [nj, ni+nj) to [begin_i, end_i) and
// [0, nj) to [begin_j, end_j) respectively.


class BipartitePairGenerator : public PairGenerator {
  const std::size_t begin_i, end_i, begin_j, end_j, ni, nj;
public:
  // Requires end_j <= begin_i when both blocks are nonempty so every
  // pair has j0 < i0. Empty generators are allowed.
  BipartitePairGenerator(
    const std::size_t begin_i,
    const std::size_t end_i,
    const std::size_t begin_j,
    const std::size_t end_j
  );

  // Half-open global index ranges for the two blocks.
  std::size_t seq_begin_i() const { return begin_i; }
  std::size_t seq_end_i() const { return end_i; }
  std::size_t seq_begin_j() const { return begin_j; }
  std::size_t seq_end_j() const { return end_j; }

  BipartitePairGenerator& operator++() override;

  std::size_t forward_map(const std::size_t value) const override;

  std::size_t reverse_map(const std::size_t value) const override;

  std::size_t i0() const override;
  std::size_t j0() const override;
};

#endif // OPTIMOTU_BIPARTITEPAIRGENERATOR_H_INCLUDED