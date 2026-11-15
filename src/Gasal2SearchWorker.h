// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_GASAL2SEARCHWORKER_H
#define OPTIMOTU_GASAL2SEARCHWORKER_H

#include "SearchWorker.h"

class Gasal2SearchWorker : public SearchWorker {
protected:
  int match = 0;
  int mismatch = 1;
  int gap_open = 0;
  int gap_extend = 1;
public:
  Gasal2SearchWorker(
    const std::vector<std::string> &query,
    const std::vector<std::string> &ref,
    const double threshold,
    const std::uint8_t threads,
    const int match = 0,
    const int mismatch = 1,
    const int gap_open = 0,
    const int gap_extend = 1,
    Dissimilarity dissimilarity = Dissimilarity::IDENTITY
  ) : SearchWorker(query, ref, threshold, threads, dissimilarity),
      match(match),
      mismatch(mismatch),
      gap_open(gap_open),
      gap_extend(gap_extend) {}
};

template<int verbose>
class Gasal2SearchWorkerImpl : public Gasal2SearchWorker {
public:
  using Gasal2SearchWorker::Gasal2SearchWorker;
  void operator()(std::size_t begin, std::size_t end) override;
  bool try_run_gpu() override;
};

#endif
