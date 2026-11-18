#include "optimotu.h"
#include "HammingSearchWorker.h"
#include "hamming_kernel.h"
extern "C" {
#include "defs.h"
}

HammingSearchWorker::HammingSearchWorker(
  const std::vector<std::string> &query,
  const std::vector<std::string> &ref,
  const double threshold,
  const std::uint8_t threads,
  const int min_overlap,
  const bool ignore_gaps
) : SearchWorker(query, ref, threshold, threads),
    pss(query, ref),
    min_overlap(min_overlap),
    ignore_gaps(ignore_gaps) {
  const int nq = static_cast<int>(query.size());
  const int nr = static_cast<int>(ref.size());
  if (optimotu_hamming_use_stripe(ignore_gaps)) {
    ref_idx.resize(static_cast<std::size_t>(nr));
    for (int j = 0; j < nr; ++j) {
      ref_idx[static_cast<std::size_t>(j)] = nq + j;
    }
    pss.fill_word_major(nq, nq + nr, ref_seq_major, ref_mask_major);
  }
}

template<int verbose>
void HammingSearchWorkerImpl<verbose>::operator()(std::size_t begin, std::size_t end) {
  std::size_t begin_i = (begin * query.size()) / threads;
  std::size_t end_i = (end * query.size()) / threads;
  OPTIMOTU_DEBUG(
    2,
    << "HammingSearchWorker thread " << begin
    << " entered; sequences [" << begin_i
    << ", "<< end_i << ")" << std::endl
  );

  std::size_t my_prealigned = 0, my_aligned = 0;
  const int nr = static_cast<int>(ref.size());
  const bool use_stripe = optimotu_hamming_use_stripe(ignore_gaps);
  std::vector<double> dists;
  std::vector<int> n_ok;
  const int nq = static_cast<int>(query.size());
  const int *tstart = nullptr;
  const int *tend = nullptr;
  if (use_stripe) {
    dists.resize(static_cast<std::size_t>(nr));
    n_ok.resize(static_cast<std::size_t>(nr));
    if (!ignore_gaps) {
      tstart = pss.start.data() + nq;
      tend = pss.end.data() + nq;
    }
  }
  for (std::size_t i = begin_i; i < end_i; i++) {
    if (!hits[i]) {
      hits[i] = std::make_unique<SearchHit>();
    }
    SearchHit & hit = *hits[i];
    if (use_stripe && nr > 0) {
      pss.score_query_vs_word_major(
        static_cast<int>(i),
        ref_seq_major.data(),
        ref_mask_major.data(),
        nr,
        nr,
        min_overlap,
        dists.data(),
        n_ok.data(),
        nullptr,
        tstart,
        tend
      );
    }
    for (std::size_t j = 0; j < ref.size(); j++) {
      double max_dist = (hit.best_dist < threshold) ? hit.best_dist : threshold;
      OPTIMOTU_DEBUG(
        4,
        << "thread" << begin
        << ": seqs " << j
        << " and " << i
        << " max relevant=" << max_dist
        << std::endl
      );
      ++my_prealigned;
      double d = use_stripe
        ? dists[j]
        : pss.dist(
            static_cast<int>(i),
            static_cast<int>(j) + static_cast<int>(query.size()),
            min_overlap,
            ignore_gaps
          );

      OPTIMOTU_DEBUG(
        4,
        << "Thread " << begin
        << ": distance=" << d
        << std::endl
      );
      if (d == 1.0) continue;
      ++my_aligned;
      if (d < max_dist) {
        OPTIMOTU_DEBUG(
          4,
          << "Thread " << begin
          << ": new best distance=" << d
          << std::endl
        );
        hit.best_dist = d;
        hit.best_ref.clear();
        hit.best_ref.push_back(j);
      } else if (d == hit.best_dist) {
        OPTIMOTU_DEBUG(
          4,
          << "Thread " << begin
          << ": tied for best distance=" << d
          << std::endl
        );
        hit.best_ref.push_back(j);
      }
#ifdef OPTIMOTU_R
      RcppThread::checkUserInterrupt();
#endif
    }
  }
  {
    std::lock_guard<std::mutex> lock(mutex);
    _prealigned += my_prealigned;
    _aligned += my_aligned;
  }
  OPTIMOTU_DEBUG(2, << "Exiting thread " << begin << std::endl);
}

template class HammingSearchWorkerImpl<0>;
template class HammingSearchWorkerImpl<1>;
template class HammingSearchWorkerImpl<2>;
template class HammingSearchWorkerImpl<3>;
template class HammingSearchWorkerImpl<4>;