#ifndef OPTIMOTU_ALLPAIRGENERATOR_H_INCLUDED
#define OPTIMOTU_ALLPAIRGENERATOR_H_INCLUDED

#include "PairGenerator.h"

// Generates all pairs {i, j} where begin <= i < j and + <= j < i
class AllPairGenerator : public PairGenerator {
private:
  std::size_t i, j = 0;
public:
  AllPairGenerator(std::size_t begin, std::size_t end) :
  PairGenerator(begin, end), i(begin == 0 ? 1 : begin) {};

  std::unique_ptr<PairGenerator> subset(
      const std::size_t begin,
      const std::size_t end
  ) const override;

  bool operator()(std::pair<std::size_t, std::size_t> & pair) override;
};

#endif // OPTIMOTU_ALLPAIRGENERATOR_H_INCLUDED
