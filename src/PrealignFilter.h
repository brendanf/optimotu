#ifndef OPTIMOTU_PREALIGN_FILTER_H_INCLUDED
#define OPTIMOTU_PREALIGN_FILTER_H_INCLUDED

#include <vector>
#include <string>
#include <memory>

class PrealignFilter {
protected:
  const std::vector<std::string> & seq;
public:
  PrealignFilter(const std::vector<std::string> & seq) : seq(seq) {};
  virtual bool operator()(
      const std::size_t i,
      const std::size_t j,
      const int max_ed,
      const int min_k,
      const int max_k
  ) const = 0;
  virtual std::unique_ptr<PrealignFilter> copy() const = 0;
};

#endif //OPTIMOTU_PREALIGN_FILTER_H_INCLUDED
