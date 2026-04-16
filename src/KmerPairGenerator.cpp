#include <unordered_set>
#include <algorithm>

#include "optimotu.h"
#include "KmerPairGenerator.h"
#include "kmer.h"

void KmerPairGenerator::update_match_index() {
  this->match_index.clear();
  do {
    const std::vector<std::uint16_t> &index = (*seq_kmer_index)[this->_i];
    if (index.size() == 0) continue;
    for (auto &kmer : index) {
      for (auto &match : (*kmer_seq_index)[kmer]) {
        if (match >= this->_i) break;
        auto entry = match_index.find(match);
        if (entry == match_index.end()) {
          match_index[match] = 1;
        } else {
          entry->second++;
        }
      }
    }
  } while (match_index.size() == 0 && ++(this->_i) < this->n);
  this->match_index_j = match_index.begin();
}

KmerPairGenerator::KmerPairGenerator(
  const std::shared_ptr<kmer_seq_index_t> kmer_seq_index,
  const std::shared_ptr<seq_kmer_index_t> seq_kmer_index,
  const double udist_threshold,
  const std::size_t offset
) : DivisiblePairGenerator(seq_kmer_index->size()),
  kmer_seq_index(kmer_seq_index),
  seq_kmer_index(seq_kmer_index),
  udist_threshold(udist_threshold),
  offset(offset) {
  update_match_index();
};

void initialize_kmer_index(
  const std::vector<std::string> & seq,
  kmer_seq_index_t & kmer_seq_index,
  seq_kmer_index_t & seq_kmer_index,
  const std::size_t begin_i,
  const std::size_t end_i
) {
  // initialize the kmer-to-sequence index
  kmer_seq_index.clear();
  kmer_seq_index.reserve(65536);
  for (std::size_t k = 0; k < 65536; k++) {
    kmer_seq_index.emplace_back();
  }
  // initialize the sequence-to-kmer index
  seq_kmer_index.clear();
  seq_kmer_index.reserve(end_i - begin_i);

  // fill the indices
  std::unordered_set<std::uint16_t> found_kmers;
  for (std::size_t i = begin_i; i < end_i; i++) {
    const std::string & s = seq[i];
    seq_kmer_index.emplace_back();
    if (s.size() > 7) {
      seq_kmer_index.back().reserve(s.size() - 7);
    } else {
      continue;
    }
    // kmer ending with current character
    std::uint16_t kmer = 0;
    // number of valid characters in current kmer
    std::size_t j = 0;

    // go character by character, updating the index
    for (auto c : s) {
      std::uint8_t newval = lookup(c);
      if (newval <= 3) {
        kmer = (kmer << 2) + newval;
      } else {
        // if we meet an ambiguous character, just reset to 0
        j = 0;
        kmer <<= 2;
      }
      if (j >= 7) {
        if (found_kmers.insert(kmer).second) {
          seq_kmer_index.back().push_back(kmer);
          kmer_seq_index[kmer].push_back(i - begin_i);
        }
      }
      ++j;
    }
    std::sort(seq_kmer_index.back().begin(), seq_kmer_index.back().end());
  }
}

KmerPairGenerator::KmerPairGenerator(
  const std::vector<std::string> & seq,
  const double udist_threshold,
  const std::size_t offset
) : DivisiblePairGenerator(seq.size()),
  kmer_seq_index(std::make_shared<kmer_seq_index_t>()),
  seq_kmer_index(std::make_shared<seq_kmer_index_t>()),
  udist_threshold(udist_threshold),
  offset(offset) {
  initialize_kmer_index(seq, *kmer_seq_index, *seq_kmer_index, offset, seq.size());
}

KmerPairGenerator& KmerPairGenerator::operator++() {
  while (true) {
    if (++(this->match_index_j) == this->match_index.end()) {
      ++(this->_i);
      update_match_index();
    }
    if (this->_i >= this->n) {
      this->has_more = false;
      return *this;
    }
    double udist = 1.0 - (double)this->match_index_j->second / (double) std::max(
      (*this->seq_kmer_index)[this->match_index_j->first].size(),
      (*this->seq_kmer_index)[this->_i].size()
    );
    if (udist > this->udist_threshold) continue;
    this->_j = this->match_index_j->first;
    return *this;
  }
}

std::size_t KmerPairGenerator::forward_map(const std::size_t value) const {
  return value + this->offset;
}
std::size_t KmerPairGenerator::reverse_map(const std::size_t value) const {
  return value - this->offset;
}
std::size_t KmerPairGenerator::i0() const {
  return this->_i + this->offset;
}
std::size_t KmerPairGenerator::j0() const {
  return this->_j + this->offset;
}

std::vector<std::unique_ptr<PairGenerator>> KmerPairGenerator::Builder::build(int verbose) const {

  // First initialize the kmer and sequence indices for each tile.
  OPTIMOTU_VERBOSE(
    1,
    << "Building " << n_subgenerators
    << " subgenerators for KmerPairGenerator with offset " << offset
    << std::endl;
  );
  std::vector<std::shared_ptr<kmer_seq_index_t>> kmer_seq_index;
  std::vector<std::shared_ptr<seq_kmer_index_t>> seq_kmer_index;
  double range = n;
  for (std::size_t tile_i = 0; tile_i < n_tiles; tile_i++) {
    std::size_t begin_i = tile_i * range / n_tiles;
    std::size_t end_i = (tile_i + 1) * range / n_tiles;
    kmer_seq_index.push_back(std::make_shared<kmer_seq_index_t>());
    seq_kmer_index.push_back(std::make_shared<seq_kmer_index_t>());
    initialize_kmer_index(
      seq,
      *kmer_seq_index.back(),
      *seq_kmer_index.back(),
      begin_i,
      end_i
    );
  }
  std::vector<std::unique_ptr<PairGenerator>> generators;
  generators.reserve(n_subgenerators);
  OPTIMOTU_VERBOSE(
    4,
    << "  Main diagonal subgenerators:" << std::endl;
  );
  // create the subgenerators along the main diagonal
  for (std::size_t tile_i = 0; tile_i < n_tiles; tile_i++) {
    std::size_t begin_i = tile_i * range / n_tiles;
    std::size_t end_i = (tile_i + 1) * range / n_tiles;
    if (begin_i == end_i) continue;
    OPTIMOTU_VERBOSE(
      4,
      << "    Subgenerator " << tile_i
      << ": " << begin_i << " to " << end_i
      << " (size = " << end_i - begin_i
      << ", offset = " << offset + begin_i
      << ")"
      << std::endl;
    );
    generators.push_back(std::make_unique<KmerPairGenerator>(
      kmer_seq_index[tile_i],
      seq_kmer_index[tile_i],
      udist_threshold,
      offset + begin_i
    ));
  }
  OPTIMOTU_VERBOSE(
    4,
    << "  Subgenerators below the diagonal:" << std::endl;
  );
  // create the subgenerators below the diagonal
  for (std::size_t tile_i = 1; tile_i < n_tiles; tile_i++) {
    std::size_t begin_i = tile_i * range / n_tiles;
    std::size_t end_i = (tile_i + 1) * range / n_tiles;
    if (begin_i == end_i) continue;
    for (std::size_t tile_j = 0; tile_j < tile_i; tile_j++) {
      std::size_t begin_j = tile_j * range / n_tiles;
      std::size_t end_j = (tile_j + 1) * range / n_tiles;
      if (begin_j == end_j) continue;
      OPTIMOTU_VERBOSE(
        4,
        << "    Subgenerator " << tile_i << "," << tile_j
        << ": " << begin_i << " to " << end_i
        << " and " << begin_j << " to " << end_j
        << " (size = " << end_i - begin_i
        << ", " << end_j - begin_j
        << ")"
        << std::endl;
      );
      generators.push_back(std::make_unique<BipartiteKmerPairGenerator>(
        begin_i, end_i, begin_j, end_j,
        udist_threshold,
        seq_kmer_index[tile_i],
        seq_kmer_index[tile_j],
        kmer_seq_index[tile_j]
      ));
    }
  }
  return generators;
}