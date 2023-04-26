#include "pad_strings.h"

#include <utility>
#include <cstring>

std::shared_ptr<char[]> pad_strings(
    const std::vector<std::string> &seq,
    std::size_t & seq_width
) {
  seq_width = 0;
  for (const auto &s : seq) {
    if (s.size() > seq_width) seq_width = s.size();
  }
  size_t lenout = seq_width * seq.size();
  char* out = new char[seq_width * lenout];
  std::fill_n(out, lenout, 'X');
  auto s = seq.begin();
  for (char* chari = out; chari < out + lenout; chari += seq_width) {
    memcpy(chari, s->c_str(), s->size());
    ++s;
  }
  return std::shared_ptr<char[]>(out);
}
