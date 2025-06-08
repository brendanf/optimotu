//SPDX-CopyrightText: (c) 2025 Brendan Furneaux
//SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_SEARCHHITS_H
#define OPTIMOTU_SEARCHHITS_H

#include <vector>
#include <string>

struct SearchHit {
  double best_dist = 1.0;
  std::vector<int> best_ref;
  virtual ~SearchHit() = default;
  virtual std::vector<std::string> & best_cigar() {
    OPTIMOTU_STOP("best_cigar() not implemented for SearchHit");
  }
};

struct SearchCigarHit : public SearchHit {
  std::vector<std::string> _best_cigar;
  virtual ~SearchCigarHit() = default;
  std::vector<std::string> & best_cigar() override {
    return _best_cigar;
  }
};

typedef std::vector<std::unique_ptr<SearchHit>> SearchHits;

#endif
