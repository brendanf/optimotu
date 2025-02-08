// SPDX-FileCopyrightText: 2025 Brendan Furneaux <brendan.furneaux@gmail.com>
// SPDX-License-Identifier: MIT

#include "AllPairGenerator.h"

std::unique_ptr<PairGenerator> AllPairGenerator::subset(
    const std::size_t begin,
    const std::size_t end
) const {
  return std::make_unique<AllPairGenerator>(begin, end);
};

bool AllPairGenerator::operator()(std::pair<std::size_t, std::size_t> & pair) {
  if (++j >= i) {
    ++i;
    j = 0;
  }
  if (i >= end) return false;
  pair.first = i;
  pair.second = j;
  return true;
}
