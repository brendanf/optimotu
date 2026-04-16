// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_ALLPAIRGENERATOR_H_INCLUDED
#define OPTIMOTU_ALLPAIRGENERATOR_H_INCLUDED

#include "PairGenerator.h"
#include "BipartitePairGenerator.h"

// Generates all pairs {i, j} where begin <= i < end and 0 <= j < i
class AllPairGenerator : public DivisiblePairGenerator {
protected:
  const std::size_t offset = 0;
public:  // constructor for use in divide() which begins with i > 1
  AllPairGenerator(const std::size_t n, const std::size_t offset = 0) :
    DivisiblePairGenerator(n), offset(offset) {};

  std::size_t forward_map(const std::size_t value) const override;

  std::size_t reverse_map(const std::size_t value) const override;

  std::size_t i0() const override;
  std::size_t j0() const override;

  AllPairGenerator& operator++() override;

  class Builder : public DivisiblePairGenerator::Builder {
  protected:
    std::size_t offset = 0;
  public:
    Builder(const std::size_t n, const std::size_t min_n_subgenerators = 1) :
      DivisiblePairGenerator::Builder(n, min_n_subgenerators) {};
    
    void set_offset(const std::size_t offset) { this->offset = offset; }

    std::vector<std::unique_ptr<PairGenerator>> build(int verbose = 0) const override;
  };
};

#endif // OPTIMOTU_ALLPAIRGENERATOR_H_INCLUDED
