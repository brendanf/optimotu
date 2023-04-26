#ifndef OPTIMOTU_PAIRGENERATOR_H_INCLUDED
#define OPTIMOTU_PAIRGENERATOR_H_INCLUDED

#include <memory>

// generates pairs {i, j} where begin <= i < end and 0 <= j < i
// implementing classes may or may not generate all such pairs
class PairGenerator {
protected:
  const std::size_t begin;
  const std::size_t end;
public:
  PairGenerator(const std::size_t begin, const std::size_t end) :
  begin(begin), end(end) {};

  // create a pair generator which works like this one, but on a different
  // range of i
  virtual std::unique_ptr<PairGenerator> subset(
      const std::size_t begin,
      const std::size_t end
  ) const = 0;

  // assign the next pair
  // return true if a pair was assigned, or false if there are no more pairs
  virtual bool operator()(std::pair<std::size_t, std::size_t> & pair) = 0;
};

#endif //OPTIMOTU_PAIRGENERATOR_H_INCLUDED
