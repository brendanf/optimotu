#include "PackedSequenceSet.h"

extern "C" {
#include "defs.h"
}

PackedSequenceSet::PackedSequenceSet(const std::vector<std::string> &seq) {
  num_seqs = seq.size();
  if (num_seqs == 0) return;
  packed_seq.reserve(seq.size());
  mask.reserve(seq.size());
  start.reserve(seq.size());
  end.reserve(seq.size());

  for (auto s : seq) {
    // if (s.size() != alen)
    //   OPTIMOTU_STOP("PackedSequenceSet: All sequences must have the same length.\n");
    alen = s.size();
    ulen = alen / NUCLEOTIDES_IN_WORD;
    if (alen > ulen * NUCLEOTIDES_IN_WORD)
      ulen++;
    mulen = alen / NUCLEOTIDES_IN_WORD / 4;
    if (alen > mulen * NUCLEOTIDES_IN_WORD / 4)
      mulen++;
    packed_seq.emplace_back(ulen);
    mask.emplace_back(mulen);
    start.emplace_back(0);
    end.emplace_back(0);
    nucleotide2binary(s.c_str(), alen, &packed_seq.back()[0],
                      &mask.back()[0], &start.back(), &end.back());
  }
}

PackedSequenceSet::PackedSequenceSet(
  const std::vector<std::string> & seq1,
  const std::vector<std::string> & seq2
) {
  num_seqs = seq1.size() + seq2.size();
  if (num_seqs == 0) return;
  packed_seq.reserve(num_seqs);
  mask.reserve(num_seqs);
  start.reserve(num_seqs);
  end.reserve(num_seqs);

  for (auto s : seq1) {
    alen = s.size();
    ulen = alen / NUCLEOTIDES_IN_WORD;
    if (alen > ulen * NUCLEOTIDES_IN_WORD)
      ulen++;
    mulen = alen / NUCLEOTIDES_IN_WORD / 4;
    if (alen > mulen * NUCLEOTIDES_IN_WORD / 4)
      mulen++;
    packed_seq.emplace_back(ulen);
    mask.emplace_back(mulen);
    start.emplace_back(0);
    end.emplace_back(0);
    nucleotide2binary(s.c_str(), alen, &packed_seq.back()[0],
                      &mask.back()[0], &start.back(), &end.back());
  }

  for (auto s : seq2) {
    alen = s.size();
    ulen = alen / NUCLEOTIDES_IN_WORD;
    if (alen > ulen * NUCLEOTIDES_IN_WORD)
      ulen++;
    mulen = alen / NUCLEOTIDES_IN_WORD / 4;
    if (alen > mulen * NUCLEOTIDES_IN_WORD / 4)
      mulen++;
    packed_seq.emplace_back(ulen);
    mask.emplace_back(mulen);
    start.emplace_back(0);
    end.emplace_back(0);
    nucleotide2binary(s.c_str(), alen, &packed_seq.back()[0],
                      &mask.back()[0], &start.back(), &end.back());
  }

}

double PackedSequenceSet::dist(const int i, const int j, const int min_overlap,
                               const bool ignore_gap) const {
  int s, e;
  s = start[j] >= start[i] ? start[j] : start[i];
  e = end[j] <= end[i] ? end[j] : end[i];
  if (s >= e) return 1.0;
  return pdistB(&packed_seq[i][0], &mask[i][0],
                &packed_seq[j][0], &mask[j][0],
                                           s, e, min_overlap);
}
