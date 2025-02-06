//SPDX-CopyrightText: (c) 2025 Brendan Furneaux
//SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_SEARCHHITS_H
#define OPTIMOTU_SEARCHHITS_H

#include <vector>

struct SearchHit {
  double best_dist = 1.0;
  std::vector<int> best_ref;
};

typedef std::vector<SearchHit> SearchHits;

#endif
