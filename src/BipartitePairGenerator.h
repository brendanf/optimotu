#ifndef OPTIMOTU_BIPARTITEPAIRGENERATOR_H_INCLUDED
#define OPTIMOTU_BIPARTITEPAIRGENERATOR_H_INCLUDED

#include "PairGenerator.h"

// Generates all pairs {i, j} where nj <= i < ni + nj and 0 <= j < nj.
// The forward mapping function maps [nj, ni+nj) to [begin_i, end_i) and
// [0, nj) to [begin_j, end_j) respectively.


class BipartitePairGenerator : public PairGenerator {
  const std::size_t begin_i, end_i, begin_j, end_j, ni, nj;
public:
  BipartitePairGenerator(
    const std::size_t begin_i,
    const std::size_t end_i,
    const std::size_t begin_j,
    const std::size_t end_j
  ) :
    PairGenerator(end_i - begin_i + end_j - begin_j),
    begin_i(begin_i), end_i(end_i),
    begin_j(begin_j), end_j(end_j),
    ni(end_i - begin_i), nj(end_j - begin_j) {
    // Internal _i runs from nj to ni+nj-1 (rows); _j from 0 to nj-1 (cols)
    _i = nj;
    has_more = (ni > 0 && nj > 0);
  };

  BipartitePairGenerator& operator++() override;

  std::size_t forward_map(const std::size_t value) const override;

  std::size_t reverse_map(const std::size_t value) const override;

  std::size_t i0() const override;
  std::size_t j0() const override;
};

#endif // OPTIMOTU_BIPARTITEPAIRGENERATOR_H_INCLUDED