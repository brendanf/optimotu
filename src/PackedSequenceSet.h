#ifndef OPTIMOTU_PACKEDSEQUENCESET_H
#define OPTIMOTU_PACKEDSEQUENCESET_H

#include <vector>
#include <string>

// Analogous to SequenceSetB in C code
struct PackedSequenceSet {
  int num_seqs, alen, ulen, mulen;
  std::vector<std::vector<unsigned long int>> packed_seq;
  std::vector<std::vector<unsigned long int>> mask;
  std::vector<int> start;
  std::vector<int> end;
  PackedSequenceSet(const std::vector<std::string> &seq);
  PackedSequenceSet(
    const std::vector<std::string> & seq1,
    const std::vector<std::string> & seq2
  );

  double dist(const int i, const int j, const int min_overlap,
              const bool ignore_gap) const;
};

#endif
