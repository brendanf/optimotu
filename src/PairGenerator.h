#ifndef OPTIMOTU_PAIRGENERATOR_H_INCLUDED
#define OPTIMOTU_PAIRGENERATOR_H_INCLUDED

#include <memory>
#include <vector>
#include <cmath>

// Generates pairs {i, j} where 0 < i < n and 0 <= j < i
// Implementing classes are not required to generate all such pairs, but they
// must implement the ++ operator to advance to the next pair, and
// should evaluate to false when there are no more pairs to generate.
class PairGenerator {
protected:
  const std::size_t n;
  bool has_more;
  std::size_t _i = 1, _j = 0;
public:
  PairGenerator(const std::size_t n) : n(n), has_more(n > 1) {};

  // Return the current indices
  std::size_t i() const { return _i; }
  std::size_t j() const { return _j; }

  // Return the indices on the original scale, if the generator has been mapped
  // to a different range.
  virtual std::size_t i0() const { return _i; }
  virtual std::size_t j0() const { return _j; }

  explicit operator bool() const {
    return has_more;
  }

  virtual PairGenerator& operator++() = 0;

  virtual ~PairGenerator() = default;

  // return the maximum value of i or j which this generator can return
  virtual std::size_t max_value() const { return n - 1; };

  // map from the range [0, max_value()] to the range [begin, end]
  virtual std::size_t forward_map(const std::size_t value) const { return value; };

  // inverse map from the range [begin, end] to the range [0, max_value()]
  virtual std::size_t reverse_map(const std::size_t value) const { return value; };
};

class DivisiblePairGenerator : public PairGenerator {
public:
  using PairGenerator::PairGenerator; // inherit constructor

  class Builder {
    protected:
    // number of items to generate pairs for
    std::size_t n;
    // number of tiles to divide the range into
    std::size_t n_tiles;
    // number of subgenerators to create
    std::size_t n_subgenerators;

    static std::size_t calculate_n_tiles(std::size_t min_n_subgenerators) {
      // The pairs are coordinates in a lower triangular matrix.
      // We can divide this into sub-matrices along the main diagonal,
      // each of which is also lower-triangular, and submatrices below the
      // main diagonal, each of which is bipartite.
      // If we divide the range evenly into n parts, then we get n submatrices
      // along the main diagonal, and n(n-1)/2 submatrices below the diagonal; for
      // a total of (n^2 + n)/2 submatrices.
      // Then, since we are targeting n_subgenerators generators,
      // n_tiles ~= 0.5*(sqrt(1 + 8 * n_subgenerators) + 1)
      return std::ceil(0.5*(std::sqrt(1.0 + 8.0*min_n_subgenerators) + 1.0));
    }

    static std::size_t calculate_n_subgenerators(std::size_t n_tiles) {
      return n_tiles * (n_tiles + 1) / 2;
    }

    public:

    Builder(const std::size_t n, const std::size_t min_n_subgenerators = 1) :
      n(n), n_tiles(Builder::calculate_n_tiles(min_n_subgenerators)), n_subgenerators(Builder::calculate_n_subgenerators(n_tiles)) {};
    // Divide the range into approximately min_n_subgenerators parts, and
    // return a vector containing one PairGenerator for each part.
    // The sub-generators should be disjoint and cover the entire range.
    // They should also be approximately equally sized, in the sense that each
    // generator returns approximately the same number of pairs.
    virtual std::vector<std::unique_ptr<PairGenerator>> build() const = 0;

    virtual ~Builder() = default;
  };
};

#endif //OPTIMOTU_PAIRGENERATOR_H_INCLUDED
