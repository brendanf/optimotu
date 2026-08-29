#include "PackedSequenceSet.h"
#include <sstream>

extern "C" {
#include "defs.h"
}

std::size_t PackedSequenceSet::estimate_bytes(std::size_t n, int alen) {
  if (n == 0 || alen <= 0) {
    return 0;
  }
  int ulen = alen / NUCLEOTIDES_IN_WORD;
  if (alen > ulen * NUCLEOTIDES_IN_WORD) {
    ulen++;
  }
  int mulen = alen / NUCLEOTIDES_IN_WORD / 4;
  if (alen > mulen * NUCLEOTIDES_IN_WORD / 4) {
    mulen++;
  }
  return n * (
    static_cast<std::size_t>(ulen) * sizeof(uint64_t) +
    static_cast<std::size_t>(mulen) * sizeof(uint64_t) +
    2 * sizeof(int)
  );
}

PackedSequenceSet::PackedSequenceSet(const SequenceSet &seq) {
  num_seqs = static_cast<int>(seq.size());
  if (num_seqs == 0) return;
  packed_seq.reserve(seq.size());
  mask.reserve(seq.size());
  start.reserve(seq.size());
  end.reserve(seq.size());
  alen = -1;
  for (int i = 0; i < num_seqs; i++) {
    const SequenceView &s = seq[static_cast<std::size_t>(i)];
    int alen_s = 0;
    for (std::size_t k = 0; k < s.size(); ++k) {
      char c = s[k];
      if ((c < 'a') || (c > 'z')) {
        alen_s++;
      }
    }
    if (alen == -1) {
      alen = alen_s;
    } else if (alen_s != alen) {
      std::stringstream ss;
      ss << "PackedSequenceSet:"
         << " All sequences must have the same length " << alen << ", but sequence "
         << i + 1 << " has length " << alen_s;
      OPTIMOTU_STOP(ss.str());
    }
    ulen = alen / NUCLEOTIDES_IN_WORD;
    if (alen > ulen * NUCLEOTIDES_IN_WORD)
      ulen++;
    mulen = alen / NUCLEOTIDES_IN_WORD / 4;
    if (alen > mulen * NUCLEOTIDES_IN_WORD / 4)
      mulen++;
    packed_seq.emplace_back(ulen, 0LL);
    mask.emplace_back(mulen, 0LL);
    start.emplace_back(0);
    end.emplace_back(0);
    nucleotide2binary(s.c_str(), alen, packed_seq.back().data(),
                      mask.back().data(), &(start.back()), &(end.back()));
  }
  verify(seq);
}

bool PackedSequenceSet::verify(const SequenceSet & seq) const {
  bool valid = true;
  if (seq.size() != static_cast<std::size_t>(num_seqs)) {
    OPTIMOTU_CERR << "mismatch between length of sequence set and number of sequences: "
                  << seq.size() << " != " << num_seqs << std::endl;
    return false;
  }
  for (int i = 0; i < num_seqs; i++) {
    const SequenceView &s = seq[static_cast<std::size_t>(i)];
    int alen_s = 0;
    for (std::size_t k = 0; k < s.size(); ++k) {
      char c = s[k];
      if ((c < 'a') || (c > 'z')) {
        alen_s++;
      }
    }
    if (alen_s != alen) return false;
    int ulen_s = alen_s / NUCLEOTIDES_IN_WORD;
    if (alen_s > ulen_s * NUCLEOTIDES_IN_WORD)
      ulen_s++;
    int mulen_s = alen_s / NUCLEOTIDES_IN_WORD / 4;
    if (alen_s > mulen_s * NUCLEOTIDES_IN_WORD / 4)
      mulen_s++;
    if (ulen_s != static_cast<int>(packed_seq[i].size())) {
      OPTIMOTU_CERR << "mismatch between length of sequence " << i
                    << " and packed bits: " << ulen_s
                    << " != " << packed_seq[i].size() << std::endl;
      valid = false;
    } else {
      int j = 0, k = 0;
      uint64_t packed_seq_i = packed_seq[i][k];
      for (std::size_t pos = 0; pos < s.size(); ++pos) {
        char c = s[pos];
        if ((c >= 'a') && (c <= 'z')) {
          continue;
        }
        if (j == 0) {
          if (k >= static_cast<int>(packed_seq[i].size())) {
            OPTIMOTU_CERR << "mismatch between length of sequence " << i
                          << " and packed bits: " << k
                          << " != " << packed_seq[i].size() << std::endl;
            valid = false;
            break;
          }
          packed_seq_i = packed_seq[i][k];
        }
        uint8_t code = (packed_seq_i >> ((NUCLEOTIDES_IN_WORD - j - 1)* 4)) & 0xfLL;
        if (code == 0) {
          if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            OPTIMOTU_CERR << "mismatch between sequence " << i
                          << " and packed bits: " << c
                          << " != " << std::hex << (uint16_t)code
                          << std::dec << std::endl;
            valid = false;
          }
        } else if (code == 0x1) {
          if (c != 'A') {
            OPTIMOTU_CERR << "mismatch between sequence " << i
                          << " and packed bits: " << c
                          << " != " << std::hex << (uint16_t)code
                          << std::dec << std::endl;
            valid = false;
          }
        } else if (code == 0x2) {
          if (c != 'C') {
            OPTIMOTU_CERR << "mismatch between sequence " << i
                          << " and packed bits: " << c
                          << " != " << std::hex << (uint16_t)code
                          << std::dec << std::endl;
            valid = false;
          }
        } else if (code == 0x4) {
          if (c != 'G') {
            OPTIMOTU_CERR << "mismatch between sequence " << i
                          << " and packed bits: " << c
                          << " != " << std::hex << (uint16_t)code
                          << std::dec << std::endl;
            valid = false;
          }
        } else if (code == 0x8) {
          if (c != 'T') {
            OPTIMOTU_CERR << "mismatch between sequence " << i
                          << " and packed bits: " << c
                          << " != " << std::hex << (uint16_t)code
                          << std::dec << std::endl;
            valid = false;
          }
        } else {
          OPTIMOTU_CERR << "invalid code in sequence " << i
                        << " at position " << k * NUCLEOTIDES_IN_WORD + j
                        << ": " << std::hex << (uint16_t)code
                        << std::dec << std::endl;
          valid = false;
        }
        j++;
        if (j >= NUCLEOTIDES_IN_WORD) {
          j = 0;
          k++;
        }
      }
    }
    if (mulen_s != static_cast<int>(mask[i].size())) {
      OPTIMOTU_CERR << "mismatch between length of sequence " << i
                    << " and mask: " << mulen_s
                    << " != " << mask[i].size() << std::endl;
      valid = false;
    } else {
      int j = 0, k = 0;
      uint64_t mask_i = mask[i][k];
      for (std::size_t pos = 0; pos < s.size(); ++pos) {
        char c = s[pos];
        if ((c >= 'a') && (c <= 'z')) {
          continue;
        }
        if (j == 0) {
          if (k >= static_cast<int>(mask[i].size())) {
            OPTIMOTU_CERR << "mismatch between length of sequence " << i
                          << " and mask: " << k
                          << " >= " << mask[i].size() << std::endl;
            valid = false;
            break;
          }
          mask_i = mask[i][k];
        }

        uint64_t mask_bit = (0x1LL << (NUCLEOTIDES_IN_WORD * 4 - j - 1)) & mask_i;
        if (mask_bit) {
          if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            OPTIMOTU_CERR << "mismatch between sequence " << i
                          << " and mask: " << c
                          << " != " << std::hex << mask_i
                          << std::dec << std::endl;
            valid = false;
          }
        } else {
          if (c == 'A' || c == 'C' || c == 'G' || c == 'T') {
            OPTIMOTU_CERR << "mismatch between sequence " << i
                          << " and mask: " << c
                          << " != " << std::hex << mask_bit
                          << std::dec << std::endl;
            valid = false;
          }
        }
        j++;
        if (j >= NUCLEOTIDES_IN_WORD * 4) {
          j = 0;
          k++;
        }
      }
    }
  }
  return valid;
}

PackedSequenceSet::PackedSequenceSet(
  const SequenceSet & seq1,
  const SequenceSet & seq2
) {
  num_seqs = static_cast<int>(seq1.size() + seq2.size());
  if (num_seqs == 0) return;
  packed_seq.reserve(static_cast<std::size_t>(num_seqs));
  mask.reserve(static_cast<std::size_t>(num_seqs));
  start.reserve(static_cast<std::size_t>(num_seqs));
  end.reserve(static_cast<std::size_t>(num_seqs));

  for (const auto &s : seq1) {
    alen = static_cast<int>(s.size());
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

  for (const auto &s : seq2) {
    alen = static_cast<int>(s.size());
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

PackedSequenceSet::PackedSequenceSet(const std::vector<std::string> &seq)
  : PackedSequenceSet(sequence_views_from_strings(seq)) {}

PackedSequenceSet::PackedSequenceSet(
  const std::vector<std::string> & seq1,
  const std::vector<std::string> & seq2
) : PackedSequenceSet(
      sequence_views_from_strings(seq1),
      sequence_views_from_strings(seq2)
    ) {}

double PackedSequenceSet::dist(const int i, const int j, const int min_overlap,
                               const bool ignore_gap) const {
  int s, e;
  s = start[j] >= start[i] ? start[j] : start[i];
  e = end[j] <= end[i] ? end[j] : end[i];
  if (s >= e) return 1.0;
  int num_ok, num_matches;
  double dist;
  if (ignore_gap) {
    pdistB(packed_seq[i].data(), mask[i].data(),
           packed_seq[j].data(), mask[j].data(),
           s, e, min_overlap, &num_ok, &num_matches, &dist);
  } else {
    pdistB2(packed_seq[i].data(), mask[i].data(),
            packed_seq[j].data(), mask[j].data(),
            s, e, min_overlap, &num_ok, &num_matches, &dist);
  }
  return dist;
}

std::pair<int, double> PackedSequenceSet::score_and_dist(
  const int i, const int j, const int min_overlap, const bool ignore_gap) const {
  int s, e;
  s = start[j] >= start[i] ? start[j] : start[i];
  e = end[j] <= end[i] ? end[j] : end[i];
  if (s >= e) return std::make_pair(0, 1.0);
  int num_ok, num_matches;
  double dist;
  if (ignore_gap) {
    pdistB(packed_seq[i].data(), mask[i].data(),
           packed_seq[j].data(), mask[j].data(),
           s, e, min_overlap, &num_ok, &num_matches, &dist);
  } else {
    pdistB2(packed_seq[i].data(), mask[i].data(),
            packed_seq[j].data(), mask[j].data(),
            s, e, min_overlap, &num_ok, &num_matches, &dist);
  }
  return std::make_pair(num_matches, dist);
}

std::tuple<bool, int, double> PackedSequenceSet::success_score_and_dist(
  const int i, const int j, const int min_overlap, const bool ignore_gap) const {
  int s, e;
  s = start[j] >= start[i] ? start[j] : start[i];
  e = end[j] <= end[i] ? end[j] : end[i];
  if (s >= e) return std::make_tuple(false, 0, 1.0);
  int result, num_ok, num_matches;
  double dist;
  if (ignore_gap) {
    result = pdistB(packed_seq[i].data(), mask[i].data(),
                    packed_seq[j].data(), mask[j].data(),
                    s, e, min_overlap, &num_ok, &num_matches, &dist);
  } else {
    result = pdistB2(packed_seq[i].data(), mask[i].data(),
                     packed_seq[j].data(), mask[j].data(),
                     s, e, min_overlap, &num_ok, &num_matches, &dist);
  }
  return std::make_tuple(result == 0, num_matches, dist);
}

std::tuple<bool, int, double, int, int, int, int> PackedSequenceSet::success_score_and_dist_gap(
  const int i, const int j, const int min_overlap, const bool ignore_gap) const {
  int s, e;
  s = start[j] >= start[i] ? start[j] : start[i];
  e = end[j] <= end[i] ? end[j] : end[i];
  if (s >= e) return std::make_tuple(false, 0, 1.0, 0, 0, 0, 0);
  int result, num_ok, num_matches, num_ins, num_del, max_ins, max_del;
  double dist;
  if (ignore_gap) {
    result = pdistB_gap(packed_seq[i].data(), mask[i].data(),
                       packed_seq[j].data(), mask[j].data(),
                       s, e, min_overlap, &num_ok, &num_matches, &dist, &num_ins, &num_del, &max_ins, &max_del);
  } else {
    result = pdistB2_gap(packed_seq[i].data(), mask[i].data(),
                       packed_seq[j].data(), mask[j].data(),
                       s, e, min_overlap, &num_ok, &num_matches, &dist, &num_ins, &num_del, &max_ins, &max_del);
  }
  return std::make_tuple(result == 0, num_matches, dist, num_ins, num_del, max_ins, max_del);
}
