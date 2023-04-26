#ifndef OPTIMOTU_SNEAKYSNAKEFILTER_H_INCLUDED
#define OPTIMOTU_SNEAKYSNAKEFILTER_H_INCLUDED

#include <memory>

extern "C" {
#include <SneakySnake.h>
}
#include "PrealignFilter.h"
#include "pad_strings.h"

class SneakySnakeFilter : public PrealignFilter {
private:
  const std::shared_ptr<char[]> pseq;
  std::size_t seq_width = 0;
  SneakySnakeFilter(const SneakySnakeFilter &obj);
public:
  SneakySnakeFilter(const std::vector<std::string> & seq);
  bool operator()(
      const std::size_t i,
      const std::size_t j,
      const int max_ed,
      const int min_k,
      const int max_k
  ) override;
  std::unique_ptr<PrealignFilter> copy() const override;
};

#endif // OPTIMOTU_SNEAKYSNAKEFILTER_H_INCLUDED
