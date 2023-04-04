#include <Rcpp.h>
#include <RcppParallel.h>
#include <RcppThread.h>
#include <bindings/cpp/WFAligner.hpp>

double distance(const std::string &a, const std::string &b, wfa::WFAligner &aligner) {
  auto status = aligner.alignEnd2End(a, b);
  if (status != wfa::WFAligner::StatusSuccessful) return 1.0;
  auto cigar = aligner.getAlignmentCigar();
  uint16_t match = 0, length = 0;
  for (char c : cigar) {
    if (c == 'M') {
      match++;
    }
    length++;
  }
  return 1.0 - double(match) / double(length);
}

//' @export
// [[Rcpp::export]]
double align(const std::string a, const std::string b,
             int match = 0, int mismatch = 1,
             int gap = 1, int extend = 0,
             int gap2 = 0, int extend2 = 0) {
  if (gap2 != 0 || extend2 != 0) {
    wfa::WFAlignerGapAffine2Pieces aligner(
        match, mismatch,
        gap, extend,
        gap2, extend2,
        wfa::WFAligner::Alignment);
    return distance(a, b, aligner);
  } else if (extend != 0) {
    wfa::WFAlignerGapAffine aligner(
        match, mismatch,
        gap, extend,
        wfa::WFAligner::Alignment);
    return distance(a, b, aligner);
  } else if (match == 0 && gap == mismatch) {
    wfa::WFAlignerEdit aligner(wfa::WFAligner::Alignment);
    return distance(a, b, aligner);
  } else if (mismatch == 0 && match == 0) {
    wfa::WFAlignerIndel aligner(wfa::WFAligner::Alignment);
    return distance(a, b, aligner);
  } else {
    wfa::WFAlignerGapLinear aligner(
        match, mismatch, gap,
        wfa::WFAligner::Alignment);
    return distance(a, b, aligner);
  }
}

// for a given kmer (implicit),
// what sequence was it found in (i),
// how many times (n),
// and at what location is the first occurrence (x)?
struct  KmerCount{
  size_t i;
  uint16_t n;
  uint16_t x;
  KmerCount& operator++() {
    this->n++;
    return *this;
  };
  KmerCount(size_t i, uint16_t n, uint16_t x) : i(i), n(n), x(x) {};
};

// for a given sequence (implicit)
// what kmer was found (kmer)
// how many times (n)
// at what location in the sequence is the first occurrence (x)
struct KmerHit{
  uint16_t kmer;
  uint16_t n;
  uint16_t x;
  KmerHit(uint16_t kmer, uint16_t n, uint16_t x) : kmer(kmer), n(n), x(x) {};
  KmerHit& operator++() {
    ++(this->n);
    return *this;
  };
  friend bool operator<(const KmerHit &a, const KmerHit &b) {
    return a.x < b.x;
  }
};

static const uint32_t debrujin32 = 0x077CB531;
static constexpr const uint8_t index32[] = {
  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
};

// ok, maybe this was too much optimization...
struct KmerBitField {
private :
  std::array<uint64_t, 1024> bits{0};
  std::array<uint32_t, 32> dirty{0};
  size_t n = 0;
public :
  void clear() {
    uint16_t bitsi = 0;
    for (uint8_t i = 0; i < 32; ++i) {
      uint32_t &j = dirty[i];
      while (j) {
        uint32_t l = (j & -j);
        // Rcpp::Rcout << "clearing field "
        //             << std::setfill('0') << std::setw(4) << bitsi + index32[((l * debrujin32) >> 27)]
        //             << " (j=" << std::hex << std::setw(8) << j
        //             << ", l=" << std::setw(8) << l << ")" << std::dec << std::endl;
        bits[bitsi + index32[((l * debrujin32) >> 27)]] = 0;
        j -= l;
      }
      dirty[i] = 0;
      bitsi+=32;
    }
    n = 0;
  };

  void insert(uint16_t k) {
    // Rcpp::Rcout << "setting " << std::hex << std::setw(4) << std::setfill('0') << k;
    uint16_t i = k >> 6;
    uint64_t j = (uint64_t)1 << (k & 0x3F);
    // Rcpp::Rcout << " as bit " << std::dec << (k & 0x3F) << " of field " << i;
    if ((bits[i] & j) == 0)  {
      // Rcpp::Rcout << " (not previously set)" << std::endl;
      ++n;
    }// else {
      // Rcpp::Rcout << " (already set)" << std::endl;
    // }
    bits[i] |= j;
    uint32_t l = (uint32_t)1 << (i & 0x1F);
    dirty[i >> 5] |= l;
  };

  bool check(uint16_t k) {
    // Rcpp::Rcout << "checking " << std::hex << std::setw(4) << k << " : " << ((bits[k >> 6] & ((uint64_t)1 << (k & 0x3F))) > 0) << std::dec << std::endl;
    return (bits[k >> 6] & ((uint64_t)1 << (k & 0x3F))) > 0;
  };

  size_t size() {
    return n;
  }
};

bool comp_by_kmer(const KmerHit &a, const KmerHit &b) {
  return a.kmer < b.kmer;
}

struct SparseDistanceMatrix {
  std::vector<size_t> &seq1, &seq2;
  std::vector<double> &dist1, &dist2;
  tthread::mutex mutex;

  SparseDistanceMatrix(
    std::vector<size_t> &seq1,
    std::vector<size_t> &seq2,
    std::vector<double> &dist1,
    std::vector<double> &dist2
  ) : seq1(seq1), seq2(seq2), dist1(dist1), dist2(dist2), mutex() {};

  void push_back(size_t s1, size_t s2, double d1, double d2) {
    seq1.push_back(s1);
    seq2.push_back(s2);
    dist1.push_back(d1);
    dist2.push_back(d2);
  }

  void append(
      std::vector<size_t> &seq1new,
      std::vector<size_t> &seq2new,
      std::vector<double> &dist1new,
      std::vector<double> &dist2new
  ) {
    mutex.lock();
    std::copy(seq1new.begin(), seq1new.end(), std::back_inserter(seq1));
    seq1new.clear();
    std::copy(seq2new.begin(), seq2new.end(), std::back_inserter(seq2));
    seq2new.clear();
    std::copy(dist1new.begin(), dist1new.end(), std::back_inserter(dist1));
    dist1new.clear();
    std::copy(dist2new.begin(), dist2new.end(), std::back_inserter(dist2));
    dist2new.clear();
    mutex.unlock();
  }
};

struct KmerAlignWorker : public RcppParallel::Worker {
  const std::vector<std::string> &seq;
  const std::vector<std::vector<KmerCount>> &kmer_seq_index;
  const std::vector<std::vector<KmerHit>> &seq_kmer_index;
  const std::vector<std::vector<KmerHit>> &seq_x_index;
  const int match, mismatch, gap, extend, gap2, extend2;
  const double dist_threshold;
  const double udist_threshold;
  const uint8_t threads;
  SparseDistanceMatrix &sdm;
  size_t &aligned;
  size_t &considered;

  KmerAlignWorker(
    const std::vector<std::string> &seq,
    const std::vector<std::vector<KmerCount>> &kmer_seq_index,
    const std::vector<std::vector<KmerHit>> &seq_kmer_index,
    const std::vector<std::vector<KmerHit>> &seq_x_index,
    const int match,
    const int mismatch,
    const int gap,
    const int extend,
    const int gap2,
    const int extend2,
    const double dist_threshold,
    const double udist_threshold,
    const uint8_t threads,
    SparseDistanceMatrix &sdm,
    size_t &aligned,
    size_t &considered
  ) : seq(seq), kmer_seq_index(kmer_seq_index), seq_kmer_index(seq_kmer_index),
  seq_x_index(seq_x_index),
  match(match), mismatch(mismatch), gap(gap), extend(extend), gap2(gap2), extend2(extend2),
  dist_threshold(dist_threshold), udist_threshold(udist_threshold), threads(threads),
  sdm(sdm), aligned(aligned), considered(considered) {};

  void operator()(std::size_t begin, std::size_t end) {
    double n = seq.size();
    double m = (n*n - 3.0*n + 2.0)/2.0;
    size_t my_considered = 0;
    size_t my_aligned = 0;
    size_t begin_i;
    std::vector<size_t> my_seq1;
    my_seq1.reserve(1000);
    std::vector<size_t> my_seq2;
    my_seq2.reserve(1000);
    std::vector<double> my_dist1;
    my_dist1.reserve(1000);
    std::vector<double> my_dist2;
    my_dist2.reserve(1000);
    KmerBitField match_hits;
    wfa::WFAlignerChoose aligner{match, mismatch, gap, extend, gap2, extend2,
                                 wfa::WFAligner::Alignment};
    if (begin == 0) {
      begin_i = 1;
    } else {
      begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
    }
    size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
    sdm.mutex.lock();
    std::cout << "Thread " << begin << " entered; sequences [" <<
      begin_i << ", "<< end_i << ")" << std::endl;
    sdm.mutex.unlock();
    for (size_t i = begin_i; i < end_i; i++) {
      const std::vector<KmerHit> &index = seq_kmer_index[i];
      if (index.size() == 0) continue;
      std::unordered_map<size_t, uint16_t> match_index;
      for (auto &kmer : index) {
        for (auto &match : kmer_seq_index[kmer.kmer]) {
          if (match.i >= i) break;
          size_t match_size;
          if (kmer.n > match.n) {
            match_size = match.n;
          } else {
            match_size = kmer.n;
          }
          auto entry = match_index.find(match.i);
          if (entry == match_index.end()) {
            match_index[match.i] = match_size;
          } else {
            entry->second += match_size;
          }
        }
      }
      if (match_index.size() == 0) continue;
      for (auto const & match : match_index) {
        ++my_considered;
        // Rcpp::Rcout << "#### seq " << i << " and " << match.first << " ####" << std::endl;
        double d1;
        if (seq[i].size() >= seq[match.first].size()) {
          d1 = 1.0 - match.second/(seq[match.first].size() - 7.0);
        } else {
          d1 = 1.0 - match.second/(seq[i].size() - 7.0);
        }
        if (d1 <= udist_threshold) {
          ++my_aligned;
          double d2 = distance(seq[match.first], seq[i], aligner);
          if (d2 <= dist_threshold) {
            my_seq1.push_back(match.first);
            my_seq2.push_back(i);
            my_dist1.push_back(d1);
            my_dist2.push_back(d2);

            if (my_seq1.size() == 1000) {
              sdm.append(my_seq1, my_seq2, my_dist1, my_dist2);
            }
          }
          RcppThread::checkUserInterrupt();
        } else {
          // Rcpp::Rcout << "giving up on comparison because udist is "
          //             << d1
          //             << std::endl;
        }
      }
    }
    if (my_seq1.size() > 0) {
      sdm.append(my_seq1, my_seq2, my_dist1, my_dist2);
    }
    sdm.mutex.lock();
    aligned += my_aligned;
    considered += my_considered;
    sdm.mutex.unlock();
  }
};

uint8_t lookup(char c) {
  if (c == 'A') {
    return 0;
  } else if (c == 'C') {
    return 1;
  } else if (c == 'G') {
    return 2;
  } else if (c == 'T') {
    return 3;
  } else {
    return 255;
  }
}

//' @export
// [[Rcpp::export]]
Rcpp::DataFrame distmx(std::vector<std::string> seq, double dist_threshold,
                       double udist_threshold,
                       int match = 1, int mismatch = 2, int gap_open = 10,
                       int gap_extend = 1, int gap_open2 = 0, int gap_extend2 = 0,
                       uint8_t threads = 1) {

  Rcpp::Rcout << "Indexing k-mers...";
  // index: for each kmer, which sequences is it found in, and how many times?
  std::vector<std::vector<KmerCount>> kmer_seq_index;
  for (size_t k = 0; k < 65536; k++) {
    kmer_seq_index.push_back({});
  }

  size_t i = 0;
  std::vector<std::vector<KmerHit>> seq_kmer_index, seq_x_index;
  std::unordered_map<uint16_t, uint16_t> kmer_location;
  // index: for each sequence, which kmers are found in it, and how many times?
  for (auto & s : seq) {
    if (s.size() > 7) {
      seq_kmer_index.emplace_back();
      seq_kmer_index.back().reserve(s.size() - 7);
      seq_x_index.emplace_back();
      seq_x_index.back().reserve(s.size() - 7);
    } else {
      seq_kmer_index.emplace_back();
      seq_x_index.emplace_back();
      continue;
    }
    uint16_t kmer = 0;
    size_t j = 0, x = 0;
    // go character by character, updating the index
    for (auto c : s) {
      uint8_t newval = lookup(c);
      if (newval <= 3) {
        kmer = (kmer << 2) + newval;
      } else {
        // if we meet an ambiguous character, just reset to 0
        j = 0;
        kmer <<= 2;
      }
      if (j >= 7) {
        // Rcpp::Rcout << "found kmer " << std::hex << kmer << " in seq " << std::dec << i;
        auto entry = kmer_location.find(kmer);
        if (entry == kmer_location.end()) {
          // Rcpp::Rcout << " (new)" << std::endl;
          seq_kmer_index.back().emplace_back(kmer, 1, x-7);
          seq_x_index.back().emplace_back(kmer, 1, x-7);
        } else {
          // Rcpp::Rcout << " (" << entry->second << ")" << std::endl;
          ++(seq_kmer_index.back()[entry->second]);
          ++(seq_x_index.back()[entry->second]);
        }
      }
      ++j;
      ++x;
    }
    std::sort(seq_kmer_index.back().begin(), seq_kmer_index.back().end(), comp_by_kmer);
    for (const auto & k : seq_kmer_index.back()) {
      kmer_seq_index[k.kmer].emplace_back(i, k.n, k.x);
    }
    ++i;
  }
  Rcpp::Rcout << "done." << std::endl;

  // outputs
  // use std library containers rather than R objects because they are much
  // faster to fill with push_back and are more thread-safe
  std::vector<size_t> seq1, seq2;
  std::vector<double> dist1, dist2;
  size_t aligned = 0, considered = 0;

  SparseDistanceMatrix sdm {seq1, seq2, dist1, dist2};
  KmerAlignWorker worker(seq, kmer_seq_index, seq_kmer_index, seq_x_index,
                        match, mismatch, gap_open, gap_extend, gap_open2,
                        gap_extend2,
                        dist_threshold, udist_threshold, threads,
                        sdm, aligned, considered);
  if (threads > 1) {
    RcppParallel::parallelFor(0, threads, worker, 1, threads);
  } else {
    worker(0, 1);
  }

  Rcpp::Rcout  << seq1.size() << " included / "
               << aligned << " aligned / "
               << considered << " considered." << std::endl;
  //   for (auto kmer : seq_index) {
  //     // Rcpp::Rcout << "looking for matches to " << std::hex << kmer.first << std::endl;
  //     auto & index = kmer_index[kmer.first];
  //     for (auto j : index) {
  //       // Rcpp::Rcout << "kmer " << std::hex << kmer.first <<
  //         // " shared with sequence " << std::dec << j.first;
  //       size_t match_size;
  //       if (j.n >= kmer.second) {
  //         match_size = kmer.second;
  //       } else {
  //         match_size = j.n;
  //       }
  //       // Rcpp::Rcout << " (" << match_size << "x) total:";
  //
  //       auto entry = seq_match.find(j.seq);
  //       if (entry == seq_match.end()) {
  //         seq_match.emplace(j.seq, match_size);
  //         // Rcpp::Rcout << match_size << std::endl;
  //       } else {
  //         entry->second += match_size;
  //         // Rcpp::Rcout << entry->second << std::endl;
  //       }
  //     }
  //     index.emplace_back(i, kmer.second);
  //   }
  // }
  //   // Rcpp::Rcout << "creating profile...";
  //   parasail_profile_t * profile = parasail_profile_create_stats_avx_256_16(
  //     // Rcpp::as<std::string>(seq[i]).c_str(), seq[i].size(),
  //     seq[i].c_str(), seq[i].size(),
  //     &parasail_dnafull
  //   );
  //   // Rcpp::Rcout << "done" << std::endl;
  //   if (threads == 1 || seq_match.size() < 10) {
  //     for (auto const & match : seq_match) {
  //       double d;
  //       if (seq[i].size() >= seq[match.first].size()) {
  //         d = 1.0 - match.second/(seq[match.first].size() - 7.0);
  //       } else {
  //         d = 1.0 - match.second/(seq[i].size() - 7.0);
  //       }
  //       if (d <= dist_threshold) {
  //         seq1.push_back(match.first);
  //         seq2.push_back(i);
  //         dist1.push_back(d);
  //         // Rcpp::Rcout << "aligning " << match.first << " and " << i << "...";
  //         // dist2.push_back(profile_align(profile, Rcpp::as<std::string>(seq[match.first])));
  //         dist2.push_back(profile_align(profile, seq[match.first]));
  //         // Rcpp::Rcout << "done" << std::endl;
  //       }
  //     }
  //   } else {
  //     std::vector<size_t> matches;
  //     matches.reserve(seq_match.size());
  //     for (auto const & match : seq_match) {
  //       double d;
  //       if (seq[i].size() >= seq[match.first].size()) {
  //         d = 1.0 - match.second/(seq[match.first].size() - 7.0);
  //       } else {
  //         d = 1.0 - match.second/(seq[i].size() - 7.0);
  //       }
  //       if (d <= dist_threshold) {
  //         seq1.push_back(match.first);
  //         seq2.push_back(i);
  //         dist1.push_back(d);
  //         matches.push_back(match.first);
  //       }
  //     }
  //     std::vector<double> dist(matches.size());
  //     AlignWorker worker(seq, matches, dist, profile);
  //     RcppParallel::parallelFor(0, matches.size(), worker, 10, threads);
  //     // Rcpp::Rcout << "adding " << dist.size() << " distances to dist2." << std::endl;
  //     for (auto const d: dist) {
  //       dist2.push_back(d);
  //     }
  //   }
  //   parasail_profile_free(profile);
  //   seq_match.clear();
  //   seq_index.clear();
  //   Rcpp::checkUserInterrupt();
  //   // Rcpp::Rcout << i << std::endl;
  //   ++i;
  // }

  Rcpp::DataFrame out = Rcpp::DataFrame::create(
    Rcpp::Named("seq1") = Rcpp::wrap(seq1),
    Rcpp::Named("seq2") = Rcpp::wrap(seq2),
    Rcpp::Named("dist1") = Rcpp::wrap(dist1),
    Rcpp::Named("dist2") = Rcpp::wrap(dist2)
  );
  return out;
}

struct PrealignAlignWorker : public RcppParallel::Worker {
  const std::vector<std::string> &seq;
  const int match, mismatch, gap, extend, gap2, extend2;
  const double dist_threshold, sim_threshold, sim_threshold_plus_1;
  const bool do_prealign, is_constrained, is_score_constrained;
  const uint8_t threads;
  SparseDistanceMatrix &sdm;
  size_t &prealigned, &aligned;

  PrealignAlignWorker(
    const std::vector<std::string> &seq,
    const int match,
    const int mismatch,
    const int gap,
    const int extend,
    const int gap2,
    const int extend2,
    const double dist_threshold,
    const bool do_prealign,
    const bool constrain,
    const uint8_t threads,
    SparseDistanceMatrix &sdm,
    size_t &prealigned,
    size_t &aligned
  ) : seq(seq),
  match(match), mismatch(mismatch), gap(gap), extend(extend), gap2(gap2), extend2(extend2),
  dist_threshold(dist_threshold), sim_threshold(1.0 - dist_threshold),
  sim_threshold_plus_1(1.0 + sim_threshold), do_prealign(do_prealign),
  is_constrained(constrain),
  is_score_constrained(
    constrain &&
      match == 0 &&
      mismatch == 1 &&
      gap == 0 &&
      extend == 1 &&
      gap2 == 0 &&
      (extend2 == 0 || extend2 == 1)
  ),
  threads(threads),
  sdm(sdm), prealigned(prealigned), aligned(aligned) {};

  void operator()(std::size_t begin, std::size_t end) {
    double n = seq.size();
    double m = (n*n - 3.0*n + 2.0)/2.0;
    size_t my_prealigned = 0;
    size_t my_aligned = 0;
    size_t begin_i;
    std::vector<size_t> my_seq1;
    my_seq1.reserve(1000);
    std::vector<size_t> my_seq2;
    my_seq2.reserve(1000);
    std::vector<double> my_dist1;
    my_dist1.reserve(1000);
    std::vector<double> my_dist2;
    my_dist2.reserve(1000);
    KmerBitField match_hits;
    wfa::WFAlignerEdit prealigner{wfa::WFAligner::Score};

    wfa::WFAlignerChoose aligner{match, mismatch, gap, extend, gap2, extend2,
                                 wfa::WFAligner::Alignment};
    if (begin == 0) {
      begin_i = 1;
    } else {
      begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
    }
    size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
    sdm.mutex.lock();
    std::cout << "Thread " << begin << " entered; sequences [" <<
      begin_i << ", "<< end_i << ")" << std::endl;
    sdm.mutex.unlock();
    for (size_t i = begin_i; i < end_i; i++) {
      for (size_t j = 0; j < i; j++) {
        double l1 = seq[i].size(), l2 = seq[j].size(), maxd1, maxs1;
        // // Rcpp::Rcout << "#### seq " << i << " (l1=" << l1 << ") and "
        // //             << j << " (l2=" << l2 <<")####" << std::endl;

        if (l1 / l2 < sim_threshold || l2/l1 < sim_threshold) continue;
        double maxd1 = dist_threshold * (l1 + l2) / sim_threshold_plus_1;
        int max_k = (int)ceil((l2 - l1 * sim_threshold) / sim_threshold_plus_1);
        int min_k = -(int)ceil((l1 - l2 * sim_threshold) / sim_threshold_plus_1);

        if (do_prealign) {
          prealigner.setMaxAlignmentScore((int) maxd1 + 1);
          // // std::cout << "(skipping) Setting band heuristics to " << min_k << ", " << max_k << std::endl;
          prealigner.setHeuristicBandedStatic(min_k, max_k);
          // std::cout << "Prealigning..." << std::endl;
          auto status = prealigner.alignEnd2End(seq[i], seq[j]);
          // std::cout << "Prealignment finished." << std::endl;
          ++my_prealigned;
          if (status != wfa::WFAligner::AlignmentStatus::StatusSuccessful) continue;
          // std::cout << "Prealignment successful." << std::endl;
          double d1 = prealigner.getAlignmentScore();
          if (d1 > maxd1) continue;
          d1 /= (l1 + l2)/sim_threshold_plus_1;
        }
        if (is_constrained) {
          aligner.setHeuristicBandedStatic(min_k, max_k);
          if (is_score_constrained) {
            aligner.setMaxAlignmentScore((int)maxd1 + 1);
          }
        }
        double d2 = distance(seq[j], seq[i], aligner);
        my_aligned++;
        if (d2 <= dist_threshold) {
          my_seq1.push_back(j);
          my_seq2.push_back(i);
          my_dist1.push_back(d1);
          my_dist2.push_back(d2);

          if (my_seq1.size() == 1000) {
            sdm.append(my_seq1, my_seq2, my_dist1, my_dist2);
          }
        }
        RcppThread::checkUserInterrupt();
      }
    }
    if (my_seq1.size() > 0) {
      sdm.append(my_seq1, my_seq2, my_dist1, my_dist2);
    }
    sdm.mutex.lock();
    aligned += my_aligned;
    prealigned += my_prealigned;
    sdm.mutex.unlock();
  }
};

//' @export
// [[Rcpp::export]]
 Rcpp::DataFrame distmx2(std::vector<std::string> seq, double dist_threshold,
                        int match = 1, int mismatch = 2, int gap_open = 10,
                        int gap_extend = 1, int gap_open2 = 0, int gap_extend2 = 0,
                        bool prealign = true, bool constrain = true,
                        uint8_t threads = 1) {
   size_t prealigned = 0, aligned = 0;

   std::vector<size_t> seq1, seq2;
   std::vector<double> dist1, dist2;

   SparseDistanceMatrix sdm {seq1, seq2, dist1, dist2};
   PrealignAlignWorker worker(seq,
                         match, mismatch, gap_open, gap_extend, gap_open2,
                         gap_extend2,
                         dist_threshold, prealign, constrain, threads,
                         sdm, prealigned, aligned);
   if (threads > 1) {
     RcppParallel::parallelFor(0, threads, worker, 1, threads);
   } else {
     worker(0, 1);
   }

   Rcpp::Rcout << seq1.size() << " included / "
               << aligned << " aligned / "
               << prealigned << " prealigned"
               << std::endl;

   Rcpp::DataFrame out = Rcpp::DataFrame::create(
     Rcpp::Named("seq1") = Rcpp::wrap(seq1),
     Rcpp::Named("seq2") = Rcpp::wrap(seq2),
     Rcpp::Named("dist1") = Rcpp::wrap(dist1),
     Rcpp::Named("dist2") = Rcpp::wrap(dist2)
   );
   return out;
}
