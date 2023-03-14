#include <Rcpp.h>
#include <parasail.h>
#include <parasail/matrices/dnafull.h>
#include <RcppParallel.h>
#include <RcppThread.h>

//' @export
// [[Rcpp::export]]
double align(std::string a, std::string b) {
  parasail_result * r = parasail_nw_banded(
    a.c_str(), a.size(),
    b.c_str(), b.size(),
    2, 1, 64,
    &parasail_dnafull);
  // int match = parasail_result_get_matches(r);
  // int length = parasail_result_get_length(r);
  // parasail_traceback_generic(
  //   a.c_str(), a.size(),
  //   b.c_str(), b.size(),
  //   "A:", "B:", &parasail_nuc44, r, '|', '*', '*', 60, 7, 1);
  parasail_result_free(r);
  return 1.0;// - double(match) / double(length);
}


//' @export
 // [[Rcpp::export]]
 double align2(const std::string &a, const std::string &b) {
   // const parasail_pfunction_info_t *info = parasail_lookup_pfunction_info("sg_striped_16");
   Rcpp::Rcout << "creating profile...";
   parasail_profile_t *profile = parasail_profile_create_stats_avx_256_16(
     a.c_str(), a.size(),
     &parasail_dnafull
   );
   Rcpp::Rcout << "done." << std::endl;

   Rcpp::Rcout << "aligning...";
   parasail_result *r = parasail_sg_stats_scan_profile_avx2_256_16(
     profile,
     b.c_str(), b.size(),
     2, 1);
   Rcpp::Rcout << "done." << std::endl;
   double match = parasail_result_get_matches(r);
   double length = parasail_result_get_length(r);
   // parasail_traceback_generic(
   //   a.c_str(), a.size(),
   //   b.c_str(), b.size(),
   //   "A:", "B:", &parasail_nuc44, r, '|', '*', '*', 60, 7, 1);
   parasail_result_free(r);
   parasail_profile_free(profile);
   return 1.0 - match / length;
 }

double profile_align(const parasail_profile_t *a, const std::string b) {
  parasail_result * r = parasail_sg_stats_scan_profile_avx2_256_16(
    a,
    b.c_str(), b.size(),
    2, 1
  );
  int match = parasail_result_get_matches(r);
  int length = parasail_result_get_length(r);
  // parasail_traceback_generic(
  //   a.c_str(), a.size(),
  //   b.c_str(), b.size(),
  //   "A:", "B:", &parasail_nuc44, r, '|', '*', '*', 60, 7, 1);
  parasail_result_free(r);
  return 1.0 - double(match) / double(length);
}

struct Anchor {
  uint16_t start1;
  uint16_t start2;
  uint16_t end1;
  uint16_t end2;
  uint16_t length;
  Anchor(uint16_t start1, uint16_t start2, uint16_t length):
    start1(start1), start2(start2), end1(start1 + length), end2(start2 + length), length(length) {};
  Anchor& operator++() {
    ++(this->end1);
    ++(this->end2);
    ++(this->length);
    return *this;
  };
  Anchor& operator+=(int i) {
    this->end1 += i;
    this->end2 += i;
    this->length += i;
    return *this;
  }
};

// struct AlignWorker : public RcppParallel::Worker {
//   const std::vector<std::string> &seq;
//   const std::vector<size_t> &matches;
//   std::vector<double> &dist;
//   const parasail_profile_t * profile;
//
//   template <class InputIt>
//   AlignWorker(
//     const std::vector<std::string> &seq,
//     const InputIt begin, const InputIt end,
//     std::vector<double> &dist,
//     const parasail_profile_t * profile
//   ): seq(seq), matches(begin, end), dist(dist), profile(profile) {};
//
//   AlignWorker(
//     const std::vector<std::string> &seq,
//     const std::vector<size_t> &matches,
//     std::vector<double> &dist,
//     const parasail_profile_t * profile
//   ): seq(seq), matches(matches), dist(dist), profile(profile) {};
//
//   void operator()(std::size_t begin, std::size_t end) {
//     for (size_t i = begin; i < end; i++) {
//       dist[i] = profile_align(profile, seq[i]);
//     }
//   };
// };

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

bool comp_by_kmer(const KmerHit &a, const KmerHit &b) {
  return a.kmer < b.kmer;
}

// check whether we can do an anchored alignment
// we will allow it if all shared kmer hits occur in the same order in each
// sequence, and if the total length of the anchor regions is 16
// (i.e., 2 different 8-mers, or a run of 9 consecutive 8mers)
bool find_anchors(const std::vector<KmerHit> &seq_kmer_index1, const std::vector<KmerHit> &seq_kmer_index2,
                  const std::vector<KmerHit> &seq_x_index1, const std::vector<KmerHit> &seq_x_index2,
                  std::vector<Anchor> &anchors) {
  auto i1 = seq_kmer_index1.begin();
  auto i2 = seq_kmer_index2.begin();
  auto end1 = seq_kmer_index1.end();
  auto end2 = seq_kmer_index2.end();
  std::unordered_set<uint16_t> match_kmers;
  while (true) {
    if (i1->kmer < i2->kmer) {
      if (++i1 == end1) break;
    } else if (i1->kmer > i2->kmer) {
      if (++i2 == end2) break;
    } else {
      if (i1->n == 1 && i2->n == 1) {
        match_kmers.insert(i1->kmer);
      }
      if (++i1 == end1) break;
      if (++i2 == end2) break;
    }
  }
  i1 = seq_x_index1.begin();
  i2 = seq_x_index2.begin();
  end1 = seq_x_index1.end();
  end2 = seq_x_index2.end();
  // total anchor length
  uint16_t l = 0, a_ext, b_ext, n = 0;
  anchors.clear();
  while(n < match_kmers.size() && i1 < end1 && i2 < end2) {
    while(i1 < end1 && match_kmers.count(i1->kmer) == 0) i1++;
    while(i2 < end2 && match_kmers.count(i2->kmer) == 0) i2++;
    if (i1->kmer == i2->kmer) {
      n++;
      if (anchors.size() == 0) {
        anchors.emplace_back(i1->x, i2->x, 8);
        // Rcpp::Rcout << "new anchor: start1=" << anchors.back().start1
        //             << " end1=" << anchors.back().end1
        //             << " start2=" << anchors.back().start2
        //             << " end2=" << anchors.back().end2
        //             << " length=" << anchors.back().length;
        l += 8;
        // Rcpp::Rcout << " total length: " << l << std::endl;
      } else if (anchors.back().end1 >= i1->x && anchors.back().end2 >= i2->x){
        a_ext = i1->x + 8 - anchors.back().end1;
        b_ext = i2->x + 8 - anchors.back().end2;
        // this happens with variable length repeats
        if (a_ext == b_ext) {
        anchors.back() += a_ext;
        l += a_ext;
        }
      } else {
        // Rcpp::Rcout << "ext anchor: start1=" << anchors.back().start1
        //             << " end1=" << anchors.back().end1
        //             << " start2=" << anchors.back().start2
        //             << " end2=" << anchors.back().end2
        //             << " length=" << anchors.back().length
        //             << " total length: " << l << std::endl;
        // anchors.emplace_back(i1->x, i2->x, 8);
        // Rcpp::Rcout << "new anchor: start1=" << anchors.back().start1
        //             << " end1=" << anchors.back().end1
        //             << " start2=" << anchors.back().start2
        //             << " end2=" << anchors.back().end2
        //             << " length=" << anchors.back().length;
        // l += 8;
        // Rcpp::Rcout << " total length: " << l << std::endl;
      }
    } else {
      // Rcpp::Rcout << "giving up on anchors because of ordering mismatch"
      //             << std::endl;
      return false;
    }
    ++i1;
    ++i2;
  }
  if (anchors.back().length > 8) {
    // Rcpp::Rcout << "ext anchor: start1=" << anchors.back().start1
    //             << " end1=" << anchors.back().end1
    //             << " start2=" << anchors.back().start2
    //             << " end2=" << anchors.back().end2
    //             << " length=" << anchors.back().length
    //             << " total length: " << l << std::endl;
  }
  if (l >= 16) return true;
  // Rcpp::Rcout << "giving up on anchors because total match is too short"
  //             << std::endl;
  return false;
}

double anchor_align(const std::string &seq1, const std::string &seq2,
                   const std::vector<Anchor> &anchors) {
  const char* s1 = seq1.c_str();
  const char* s2 = seq2.c_str();
  uint16_t x1 = 0, x2 = 0;
  double matches = 0.0, length = 0.0;
  auto a_i = anchors.begin();
  auto a_end = anchors.end();
  if (a_i->start1 == 0 || a_i->start2 == 0) {
    matches += a_i->length;
    length += a_i->length;
  } else {
    // Rcpp::Rcout << "aligning seq 1:[0," << a_i->start1 << ") with seq 2:[0," <<
    //   a_i->start2 << ")" << std::endl;
    parasail_result_t *r = parasail_sg_qb_db_stats_scan_avx2_256_16(
      s1, a_i->start1,
      s2, a_i->start2,
      2, 1,
      &parasail_dnafull
    );
    matches += parasail_result_get_matches(r) + a_i->length;
    length += parasail_result_get_length(r) + a_i->length;
    parasail_result_free(r);
  }
  x1 = a_i->end1;
  x2 = a_i->end2;
  ++a_i;
  for (;a_i != a_end; ++a_i) {
    uint32_t l1 = a_i->start1 - x1;
    uint32_t l2 = a_i->start2 - x2;
    if (l1 == 0) {
      length += l2;
      x1 = a_i->end1;
      x2 = a_i->end2;
      continue;
    }
    if (l2 == 0) {
      length += l1;
      x1 = a_i->end1;
      x2 = a_i->end2;
      continue;
    }

    // Rcpp::Rcout << "aligning seq 1:[" << x1 << "," << a_i->start1 <<
    //   ") with seq 2:[" << x2 << "," << a_i->start2 << ")" << std::endl;
    parasail_result *r = parasail_sw_stats_scan_avx2_256_16(
      s1 + x1, l1,
      s2 + x2, l2,
      2, 1,
      &parasail_dnafull
    );
    x1 = a_i->end1;
    x2 = a_i->end2;
    matches += parasail_result_get_matches(r) + a_i->length;
    length += parasail_result_get_length(r) + a_i->length;
    parasail_result_free(r);
  }
  --a_i;
  if (a_i->end1 < seq1.size() && a_i->end2 < seq2.size()) {
    parasail_result *r = parasail_sg_qe_de_stats_scan_avx2_256_16(
      s1 + x1, seq1.size() - x1,
      s2 + x2, seq2.size() - x2,
      2, 1,
      &parasail_dnafull
    );
    matches += parasail_result_get_matches(r) + a_i->length;
    length += parasail_result_get_length(r) + a_i->length;
    parasail_result_free(r);
  }
  return 1.0 - (matches/length);
}

struct BigAlignWorker : public RcppParallel::Worker {
  const std::vector<std::string> &seq;
  const std::vector<std::vector<KmerCount>> &kmer_seq_index;
  const std::vector<std::vector<KmerHit>> &seq_kmer_index;
  const std::vector<std::vector<KmerHit>> &seq_x_index;
  const double dist_threshold;
  const bool heuristic;
  const uint8_t threads;
  std::vector<size_t> &seq1, &seq2;
  std::vector<double> &dist1, &dist2;
  tthread::mutex &mutex;
  size_t &aligned;
  size_t &anchored;

  BigAlignWorker(
    const std::vector<std::string> &seq,
    const std::vector<std::vector<KmerCount>> &kmer_seq_index,
    const std::vector<std::vector<KmerHit>> &seq_kmer_index,
    const std::vector<std::vector<KmerHit>> &seq_x_index,
    const double dist_threshold,
    const uint8_t threads,
    const bool heuristic,
    std::vector<size_t> &seq1,
    std::vector<size_t> &seq2,
    std::vector<double> &dist1,
    std::vector<double> &dist2,
    tthread::mutex &mutex,
    size_t &aligned,
    size_t &anchored
  ) : seq(seq), kmer_seq_index(kmer_seq_index), seq_kmer_index(seq_kmer_index),
  seq_x_index(seq_x_index),
  dist_threshold(dist_threshold), threads(threads), heuristic(heuristic),
  seq1(seq1), seq2(seq2), dist1(dist1), dist2(dist2),
  mutex(mutex), aligned(aligned), anchored(anchored) {};

  void operator()(std::size_t begin, std::size_t end) {
    double n = seq.size();
    double m = (n*n - 3.0*n + 2.0)/2.0;
    size_t my_anchored = 0;
    size_t my_aligned = 0;
    size_t begin_i;
    std::vector<size_t> my_seq1;
    my_seq1.reserve(100);
    std::vector<size_t> my_seq2;
    my_seq2.reserve(100);
    std::vector<double> my_dist1;
    my_dist1.reserve(100);
    std::vector<double> my_dist2;
    my_dist2.reserve(100);

    if (begin == 0) {
      begin_i = 1;
    } else {
      begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
    }
    size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
    mutex.lock();
    std::cout << "Thread " << begin << " entered; sequences [" <<
      begin_i << ", "<< end_i << ")" << std::endl;
    mutex.unlock();
    std::vector<Anchor> anchors;
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
      parasail_profile_t * profile = parasail_profile_create_stats_avx_256_16(
        // Rcpp::as<std::string>(seq[i]).c_str(), seq[i].size(),
        seq[i].c_str(), seq[i].size(),
        &parasail_dnafull
      );
      for (auto const & match : match_index) {
        // Rcpp::Rcout << "#### seq " << i << " and " << match.first << " ####" << std::endl;
        double d1;
        if (seq[i].size() >= seq[match.first].size()) {
          d1 = 1.0 - match.second/(seq[match.first].size() - 7.0);
        } else {
          d1 = 1.0 - match.second/(seq[i].size() - 7.0);
        }
        if (d1 <= dist_threshold) {
          double d2;
          if (
              heuristic
              && find_anchors(seq_kmer_index[i], seq_kmer_index[match.first],
                           seq_x_index[i], seq_x_index[match.first], anchors)
          ) {
            d2 = anchor_align(seq[i], seq[match.first], anchors);
            anchors.clear();
            ++my_aligned;
            ++my_anchored;
          } else {
            // Rcpp::Rcout << "profile aligning full length sequences" << std::endl;
            d2 = profile_align(profile, seq[match.first]);
            ++my_aligned;
          }
          // double d2 = align(seq[i], seq[match.first]);

          my_seq1.push_back(match.first);
          my_seq2.push_back(i);
          my_dist1.push_back(d1);
          my_dist2.push_back(d2);
          if (my_seq1.size() == 100) {
            mutex.lock();
            std::copy(my_seq1.begin(), my_seq1.end(), std::back_inserter(seq1));
            my_seq1.clear();
            std::copy(my_seq2.begin(), my_seq2.end(), std::back_inserter(seq2));
            my_seq2.clear();
            std::copy(my_dist1.begin(), my_dist1.end(), std::back_inserter(dist1));
            my_dist1.clear();
            std::copy(my_dist2.begin(), my_dist2.end(), std::back_inserter(dist2));
            my_dist2.clear();
            mutex.unlock();
          }
          RcppThread::checkUserInterrupt();
        } else {
          // Rcpp::Rcout << "giving up on comparison because udist is "
          //             << d1
          //             << std::endl;
        }
      }
      parasail_profile_free(profile);
    }
    mutex.lock();
    aligned += my_aligned;
    anchored += my_anchored;
    mutex.unlock();
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
                       uint8_t threads = 1, bool heuristic = true) {

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
  size_t aligned = 0, anchored = 0;

  tthread::mutex mutex;
  BigAlignWorker worker(seq, kmer_seq_index, seq_kmer_index, seq_x_index,
                        dist_threshold, threads, heuristic,
                        seq1, seq2, dist1, dist2, mutex, aligned, anchored);
  RcppParallel::parallelFor(0, threads, worker, 1, threads);

  Rcpp::Rcout << anchored << "/" << aligned << " used anchored alignment." << std::endl;
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
