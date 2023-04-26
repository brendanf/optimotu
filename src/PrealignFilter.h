#ifndef OPTIMOTU_PREALIGN_FILTER_H_INCLUDED
#define OPTIMOTU_PREALIGN_FILTER_H_INCLUDED

#include <vector>
#include <string>

class PrealignFilter {
protected:
  const std::vector<std::string> & seq;
public:
  PrealignFilter(const std::vector<std::string> & seq) : seq(seq) {};
  virtual bool operator()(std::size_t i, std::size_t j, int d) = 0;
};

#endif //OPTIMOTU_PREALIGN_FILTER_H_INCLUDED
