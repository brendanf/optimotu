#include "optimotu.h"
#include "BipartiteKmerPairGenerator.h"
#include <algorithm>
#include <cmath>

void BipartiteKmerPairGenerator::update_match_index() {
  this->match_index.clear();
  do {
    const std::vector<std::uint16_t> &index = (*this->seq_kmer_index_i)[this->_i - this->nj];
    if (index.size() == 0) continue;
    for (auto &kmer : index) {
      for (auto &match : (*this->kmer_seq_index_j)[kmer]) {
        auto entry = this->match_index.find(match);
        if (entry == this->match_index.end()) {
          this->match_index[match] = 1;
        } else {
          entry->second++;
        }
      }
    }
  } while (match_index.size() == 0 && ++(this->_i) < this->ni + this->nj);
  this->match_index_j = match_index.begin();
}

BipartiteKmerPairGenerator& BipartiteKmerPairGenerator::operator++() {
  while (this->has_more) {
    if (this->match_index_j == this->match_index.end()) {
      ++(this->_i);
      if (this->_i >= this->ni + this->nj) {
        this->has_more = false;
        return *this;
      }
      update_match_index();
      if (this->match_index.empty()) continue;
      this->match_index_j = this->match_index.begin();
    } else {
      ++(this->match_index_j);
      if (this->match_index_j == this->match_index.end()) continue;
    }
    if (this->_i >= this->ni + this->nj) {
      this->has_more = false;
      return *this;
    }
    double udist = 1.0 - (double)this->match_index_j->second / (double) std::max(
      (*this->seq_kmer_index_j)[this->match_index_j->first].size(),
      (*this->seq_kmer_index_i)[this->_i - this->nj].size()
    );
    if (udist > this->udist_threshold) continue;
    this->_j = this->match_index_j->first;
    return *this;
  }
  return *this;
}

std::size_t BipartiteKmerPairGenerator::forward_map(const std::size_t value) const {
  if (value < this->nj) {
    return value + this->begin_j;
  } else if (value < this->ni + this->nj) {
    return value - this->nj + this->begin_i;
  } else {
    OPTIMOTU_STOP(
      "Invalid value: " + std::to_string(value) +
      " for BipartiteKmerPairGenerator::forward_map() with ni = " + std::to_string(ni) +
      " and nj = " + std::to_string(nj)
    );
  }
}

std::size_t BipartiteKmerPairGenerator::reverse_map(const std::size_t value) const {
  if (value >= this->begin_j && value < this->end_j) {
    return value - this->begin_j;
  } else if (value >= this->begin_i && value < this->end_i) {
    return value - this->begin_i + this->nj;
  } else {
    OPTIMOTU_STOP(
      "Invalid value: " + std::to_string(value) +
      " for BipartiteKmerPairGenerator::reverse_map() with ni = " + std::to_string(ni) +
      " and nj = " + std::to_string(nj)
    );
  }
}

std::size_t BipartiteKmerPairGenerator::i0() const {
  return this->_i - this->nj + this->begin_i;
}
std::size_t BipartiteKmerPairGenerator::j0() const {
  return this->_j + this->begin_j;
}
std::size_t BipartiteKmerPairGenerator::max_value() const {
  return this->ni + this->nj - 1;
}