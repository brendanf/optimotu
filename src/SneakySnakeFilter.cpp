#include "SneakySnakeFilter.h"

SneakySnakeFilter::SneakySnakeFilter(const SneakySnakeFilter &obj) :
  PrealignFilter(obj.seq), pseq(obj.pseq), seq_width(obj.seq_width) {};

SneakySnakeFilter::SneakySnakeFilter(const std::vector<std::string> & seq) :
  PrealignFilter(seq), pseq(pad_strings(seq, seq_width)) {};

bool SneakySnakeFilter::operator()(
    const std::size_t i,
    const std::size_t j,
    const int max_ed,
    const int min_k,
    const int max_k
) const {
  return (bool) SneakySnake(
    seq[j].size(),
    pseq.get() + j*seq_width,
    pseq.get() + i*seq_width,
    max_ed,
    2*max_k + 1,
    0,
    seq[j].size()
  );
}

std::unique_ptr<PrealignFilter> SneakySnakeFilter::copy() const {
  return std::unique_ptr<SneakySnakeFilter>(new SneakySnakeFilter(*this));
}
