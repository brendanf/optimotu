// SPDX-FileCopyrightText: 2025 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_SEQUENCEVIEW_H_INCLUDED
#define OPTIMOTU_SEQUENCEVIEW_H_INCLUDED

#include <cstddef>
#include <string>
#include <vector>

// Non-owning view of a contiguous ASCII sequence buffer. Safe for worker
// threads when the underlying storage (R CHARSXP or std::string) outlives
// the view. Does not call the R API.
struct SequenceView {
  const char *data_ = nullptr;
  std::size_t size_ = 0;

  SequenceView() = default;
  SequenceView(const char *data, std::size_t size) : data_(data), size_(size) {}
  SequenceView(const std::string &s) : data_(s.data()), size_(s.size()) {}

  const char *data() const { return data_; }
  const char *c_str() const { return data_; }
  std::size_t size() const { return size_; }
  bool empty() const { return size_ == 0; }

  char operator[](std::size_t i) const { return data_[i]; }

  const char *begin() const { return data_; }
  const char *end() const { return data_ + size_; }
};

using SequenceSet = std::vector<SequenceView>;

// Build views that alias owned std::string buffers. The strings must outlive
// the returned views.
inline SequenceSet sequence_views_from_strings(
  const std::vector<std::string> &seq
) {
  SequenceSet out;
  out.reserve(seq.size());
  for (const auto &s : seq) {
    out.emplace_back(s);
  }
  return out;
}

#endif // OPTIMOTU_SEQUENCEVIEW_H_INCLUDED
