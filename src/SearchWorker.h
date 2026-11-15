// SPDX-CopyrightText: (c) 2025 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_SEARCHWORKER_H
#define OPTIMOTU_SEARCHWORKER_H

#include <cstdint>
#include <mutex>
#include <RcppParallel.h>
#include "SearchHits.h"
#include "alignment_enums.h"

class SearchWorker : public RcppParallel::Worker{
protected:
  const std::vector<std::string> & query;
  const std::vector<std::string> & ref;
  const double threshold;
  const std::uint8_t threads;
  std::mutex mutex;
  size_t _prealigned = 0, _aligned = 0;
  SearchHits hits;
  const Dissimilarity dissimilarity;
public:
  SearchWorker(
    const std::vector<std::string> & query,
    const std::vector<std::string> & ref,
    const double threshold,
    const std::uint8_t threads,
    Dissimilarity dissimilarity = Dissimilarity::IDENTITY
  ) : query(query), ref(ref), threshold(threshold), threads(threads),
  hits(query.size()), dissimilarity(dissimilarity) {}

  size_t prealigned() const { return _prealigned; }

  size_t aligned() const { return _aligned; }

  std::uint8_t n_threads() const { return threads; }

  const SearchHits & get_hits() const { return hits; }

  virtual bool try_run_gpu() { return false; }
};

#endif
