// SPDX-FileCopyrightText: 2026 Brendan Furneaux
// SPDX-License-Identifier: MIT

#ifndef OPTIMOTU_KSW2DISTWORKER_H
#define OPTIMOTU_KSW2DISTWORKER_H

#ifdef OPTIMOTU_R
#include "DistWorker.h"
#include "alignment_enums.h"

class Ksw2DistWorker : public DistWorker {
protected:
  int match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2;
public:
  Ksw2DistWorker(
    const std::vector<std::string> &seq,
    const double dist_threshold,
    DivisiblePairGenerator::Builder & pgb,
    SparseDistanceMatrix &sdm,
    int match, int mismatch,
    int gap_open, int gap_extend,
    int gap_open2, int gap_extend2,
    int verbose = 0
  );
};

template<int verbose, bool is_constrained, enum AlignmentSpan span = AlignmentSpan::GLOBAL, typename SparseDistanceMatrixType = SparseDistanceMatrix>
class Ksw2DistWorkerImpl : public Ksw2DistWorker {
  using DistWorker::seq;
  using DistWorker::dist_threshold;
  using DistWorker::threads;
  using DistWorker::sdm;
  using DistWorker::sim_threshold;
  using Ksw2DistWorker::match;
  using Ksw2DistWorker::mismatch;
  using Ksw2DistWorker::gap_open;
  using Ksw2DistWorker::gap_extend;
  using Ksw2DistWorker::gap_open2;
  using Ksw2DistWorker::gap_extend2;
public:
  Ksw2DistWorkerImpl(
    const std::vector<std::string> &seq,
    const double dist_threshold,
    DivisiblePairGenerator::Builder & pgb,
    SparseDistanceMatrix &sdm,
    int match, int mismatch,
    int gap_open, int gap_extend,
    int gap_open2, int gap_extend2
  ) : Ksw2DistWorker(seq, dist_threshold, pgb, sdm, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, verbose) {};

  virtual void operator()(std::size_t begin, std::size_t end) override;
};

std::unique_ptr<Ksw2DistWorker> create_ksw2_dist_worker(
  const std::vector<std::string> &seq,
  const double dist_threshold,
  DivisiblePairGenerator::Builder & pgb,
  SparseDistanceMatrix &sdm,
  int match = 0,
  int mismatch = 1,
  int gap_open = 0,
  int gap_extend = 1,
  int gap_open2 = 0,
  int gap_extend2 = 1,
  int verbose = 0,
  enum AlignmentSpan span = AlignmentSpan::GLOBAL,
  bool constrain = true
);

#endif //OPTIMOTU_R

#endif //OPTIMOTU_KSW2DISTWORKER_H
