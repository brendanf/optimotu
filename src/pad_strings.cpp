#include "pad_strings.h"

std::pair<char*, size_t> pad_strings(const std::vector<std::string> &seq) {
  size_t max_len = 0;
  for (const auto &s : seq) {
    if (s.size() > max_len) max_len = s.size();
  }
  size_t lenout = max_len * seq.size();
  char* out = new char[max_len * lenout];
  std::fill_n(out, lenout, 'X');
  auto s = seq.begin();
  for (char* chari = out; chari < out + lenout; chari += max_len) {
    memcpy(chari, s->c_str(), s->size());
    ++s;
  }
  return std::pair<char*, size_t>{out, max_len};
}
