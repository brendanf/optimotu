// SPDX-CopyrightText: (c) 2025 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_SEARCHWORKER_H
#define OPTIMOTU_SEARCHWORKER_H

#include <cstdint>
#include <mutex>
#include <RcppParallel.h>
#include "SearchHits.h"

class SearchWorker : public RcppParallel::Worker{
protected:
  const std::vector<std::string> & query;
  const std::vector<std::string> & ref;
  const double threshold;
  const std::uint8_t threads;
  std::mutex mutex;
  size_t _prealigned = 0, _aligned = 0;
  SearchHits hits;
public:
  SearchWorker(
    const std::vector<std::string> & query,
    const std::vector<std::string> & ref,
    const double threshold,
    const std::uint8_t threads
  ) : query(query), ref(ref), threshold(threshold), threads(threads),
  hits(query.size()) {}

  size_t prealigned() const { return _prealigned; }

  size_t aligned() const { return _aligned; }

  std::uint8_t n_threads() const { return threads; }

  const SearchHits & get_hits() const { return hits; }
};

#endif
