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
   int match = parasail_result_get_matches(r);
   int length = parasail_result_get_length(r);
   // parasail_traceback_generic(
   //   a.c_str(), a.size(),
   //   b.c_str(), b.size(),
   //   "A:", "B:", &parasail_nuc44, r, '|', '*', '*', 60, 7, 1);
   parasail_result_free(r);
   parasail_profile_free(profile);
   return 1.0 - double(match) / double(length);
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

struct AlignWorker : public RcppParallel::Worker {
  const std::vector<std::string> &seq;
  const std::vector<size_t> &matches;
  std::vector<double> &dist;
  const parasail_profile_t * profile;

  template <class InputIt>
  AlignWorker(
    const std::vector<std::string> &seq,
    const InputIt begin, const InputIt end,
    std::vector<double> &dist,
    const parasail_profile_t * profile
  ): seq(seq), matches(begin, end), dist(dist), profile(profile) {};

  AlignWorker(
    const std::vector<std::string> &seq,
    const std::vector<size_t> &matches,
    std::vector<double> &dist,
    const parasail_profile_t * profile
  ): seq(seq), matches(matches), dist(dist), profile(profile) {};

  void operator()(std::size_t begin, std::size_t end) {
    for (size_t i = begin; i < end; i++) {
      dist[i] = profile_align(profile, seq[i]);
    }
  };
};

struct  KmerCount{
  size_t i;
  size_t n;
  KmerCount& operator++() {
    this->n++;
    return *this;
  };
  KmerCount(size_t i, size_t n) : i(i), n(n) {};
};


struct BigAlignWorker : public RcppParallel::Worker {
  double dist_threshold;
  const std::vector<std::string> &seq;
  const std::vector<std::vector<KmerCount>> &kmer_index;
  const std::vector<std::unordered_map<uint16_t, uint16_t>> &seq_index;
  uint8_t threads;
  std::vector<size_t> &seq1, &seq2;
  std::vector<double> &dist1, &dist2;
  tthread::mutex &mutex;

  BigAlignWorker(
    const std::vector<std::string> &seq,
    const std::vector<std::vector<KmerCount>> &kmer_index,
    const std::vector<std::unordered_map<uint16_t, uint16_t>> &seq_index,
    double dist_threshold,
    uint8_t threads,
    std::vector<size_t> &seq1,
    std::vector<size_t> &seq2,
    std::vector<double> &dist1,
    std::vector<double> &dist2,
    tthread::mutex &mutex
  ) : seq(seq), kmer_index(kmer_index), seq_index(seq_index),
  dist_threshold(dist_threshold), threads(threads),
  seq1(seq1), seq2(seq2), dist1(dist1), dist2(dist2),
  mutex(mutex) {};

  void operator()(std::size_t begin, std::size_t end) {
    double n = seq.size();
    double m = (n*n - 3.0*n + 2.0)/2.0;
    size_t begin_i;
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
    for (size_t i = begin_i; i < end_i; i++) {
      const std::unordered_map<uint16_t, uint16_t> &index = seq_index[i];
      if (index.size() == 0) continue;
      std::unordered_map<size_t, uint16_t> match_index;
      for (auto &kmer : index) {
        for (auto &match : kmer_index[kmer.first]) {
          if (match.i >= i) break;
          size_t match_size;
          if (kmer.second > match.n) {
            match_size = match.n;
          } else {
            match_size = kmer.second;
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
        double d1;
        if (seq[i].size() >= seq[match.first].size()) {
          d1 = 1.0 - match.second/(seq[match.first].size() - 7.0);
        } else {
          d1 = 1.0 - match.second/(seq[i].size() - 7.0);
        }
        if (d1 <= dist_threshold) {
          double d2 = profile_align(profile, seq[match.first]);
          // double d2 = align(seq[i], seq[match.first]);
          mutex.lock();
          seq1.push_back(match.first);
          seq2.push_back(i);
          dist1.push_back(d1);
          dist2.push_back(d2);
          mutex.unlock();
        }
      }
    }
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
Rcpp::DataFrame distmx(std::vector<std::string> seq, double dist_threshold, uint8_t threads = 1) {

  Rcpp::Rcout << "Indexing k-mers...";
  // index: for each kmer, which sequences is it found in, and how many times?
  std::vector<std::vector<KmerCount>> kmer_index;
  for (size_t k = 0; k < 65536; k++) {
    kmer_index.push_back({});
  }

  size_t i = 0;
  std::vector<std::unordered_map<uint16_t, uint16_t>> seq_index;
  // index: for each sequence, which kmers are found in it, and how many times?
  for (auto & s : seq) {
    if (s.size() > 7) {
      seq_index.emplace_back(s.size() - 7);
    } else {
      seq_index.emplace_back(0);
      continue;
    }
    uint16_t kmer = 0;
    size_t j = 0;
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
        auto entry = seq_index.back().find(kmer);
        if (entry == seq_index.back().end()) {
          // Rcpp::Rcout << " (new)" << std::endl;
          seq_index.back().emplace(kmer, 1);
        } else {
          // Rcpp::Rcout << " (" << entry->second << ")" << std::endl;
          entry->second++;
        }
      }
      ++j;
    }
    for (const auto & k : seq_index.back()) {
      kmer_index[k.first].emplace_back(i, k.second);
    }
    ++i;
  }
  Rcpp::Rcout << "done." << std::endl;

  // outputs
  // use std library containers rather than R objects because they are much
  // faster to fill with push_back and are more thread-safe
  std::vector<size_t> seq1, seq2;
  std::vector<double> dist1, dist2;

  tthread::mutex mutex;
  BigAlignWorker worker(seq, kmer_index, seq_index, dist_threshold, threads,
                        seq1, seq2, dist1, dist2, mutex);
  RcppParallel::parallelFor(0, threads, worker, 1, threads);

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
