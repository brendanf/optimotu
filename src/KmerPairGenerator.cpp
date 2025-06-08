#include <unordered_set>
#include <algorithm>

#include "KmerPairGenerator.h"
#include "kmer.h"

void KmerPairGenerator::update_match_index() {
  this->match_index.clear();
  do {
    const std::vector<std::uint16_t> &index = (*seq_kmer_index)[this->i];
    if (index.size() == 0) continue;
    for (auto &kmer : index) {
      for (auto &match : (*kmer_seq_index)[kmer]) {
        if (match >= this->i) break;
        auto entry = match_index.find(match);
        if (entry == match_index.end()) {
          match_index[match] = 1;
        } else {
          entry->second++;
        }
      }
    }
  } while (match_index.size() == 0 && ++(this->i) < this->end);
  this->j = match_index.begin();
}

KmerPairGenerator::KmerPairGenerator(
  const std::size_t begin,
  const std::size_t end,
  const std::shared_ptr<std::vector<std::vector<std::size_t>>> kmer_seq_index,
  const std::shared_ptr<std::vector<std::vector<std::uint16_t>>> seq_kmer_index,
  const double udist_threshold
) : PairGenerator(begin, end), kmer_seq_index(kmer_seq_index),
  seq_kmer_index(seq_kmer_index), udist_threshold(udist_threshold),
  i(begin == 0 ? 1 : begin) {
  update_match_index();
};

KmerPairGenerator::KmerPairGenerator(
  const std::size_t begin,
  const std::size_t end,
  const std::vector<std::string> & seq,
  const double udist_threshold
) : PairGenerator(begin, end),
kmer_seq_index(std::make_shared<std::vector<std::vector<std::size_t>>>()),
seq_kmer_index(std::make_shared<std::vector<std::vector<std::uint16_t>>>()),
udist_threshold(udist_threshold),
i(begin == 0 ? 1 : begin) {
  // Rcpp::Rcout << "Indexing k-mers...";
  // index: for each kmer, which sequences is it found in, and how many times?
  this->kmer_seq_index->reserve(65536);
  for (std::size_t k = 0; k < 65536; k++) {
    kmer_seq_index->emplace_back();
  }

  this->seq_kmer_index->reserve(seq.size());
  std::size_t ii = 0;
  std::unordered_set<std::uint16_t> kmer_location;
  // index: for each sequence, which kmers are found in it?
  for (auto & s : seq) {
      seq_kmer_index->emplace_back();
    if (s.size() > 7) {
      seq_kmer_index->back().reserve(s.size() - 7);
    } else {
      continue;
    }
    std::uint16_t kmer = 0;
    std::size_t jj = 0;
    // go character by character, updating the index
    for (auto c : s) {
      std::uint8_t newval = lookup(c);
      if (newval <= 3) {
        kmer = (kmer << 2) + newval;
      } else {
        // if we meet an ambiguous character, just reset to 0
        jj = 0;
        kmer <<= 2;
      }
      if (jj >= 7) {
        // Rcpp::Rcout << "found kmer " << std::hex << kmer << " in seq " << std::dec << i;
        if (kmer_location.insert(kmer).second) {
          // Rcpp::Rcout << " (new)" << std::endl;
          this->seq_kmer_index->back().push_back(kmer);
          (*this->kmer_seq_index)[kmer].push_back(ii);
        }
      }
      ++jj;
    }
    std::sort(
      this->seq_kmer_index->back().begin(),
      this->seq_kmer_index->back().end()
    );
    ++ii;
  }
}

std::unique_ptr<PairGenerator> KmerPairGenerator::subset(
    const std::size_t begin,
    const std::size_t end
) const {
  return std::unique_ptr<KmerPairGenerator>(
    new KmerPairGenerator(
        begin,
        end,
        this->kmer_seq_index,
        this->seq_kmer_index,
        this->udist_threshold
    )
  );
}

bool KmerPairGenerator::operator()(std::pair<std::size_t, std::size_t> & pair) {
  while (true) {
    if (++(this->j) == this->match_index.end()) {
      ++(this->i);
      update_match_index();
    }
    if (i >= end) return false;
    double udist = 1.0 - (double)this->j->second / (double) std::max(
      (*this->kmer_seq_index)[this->j->first].size(),
      (*this->kmer_seq_index)[this->i].size()
    );
    if (udist > this->udist_threshold) continue;
    pair.first = i;
    pair.second = this->j->first;
    return true;
  }
}
