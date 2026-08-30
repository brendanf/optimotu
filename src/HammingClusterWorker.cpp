#include "HammingClusterWorker.h"
#include "AllPairGenerator.h"
#include "BipartitePairGenerator.h"
#include "MappedClusterAlgorithm.h"
#include "MultipleClusterAlgorithm.h"
extern "C" {
#include "defs.h"
}

typedef RcppParallel::RMatrix<int> matrix_t;

namespace
{

  // Hamming visits ~1e9 pairs; checking R interrupt every pair dominates.
  constexpr std::size_t kHammingInterruptBatch = 65536;

  inline void hamming_maybe_interrupt(std::size_t &pairs_since_check)
  {
    if (++pairs_since_check >= kHammingInterruptBatch)
    {
      pairs_since_check = 0;
#ifdef OPTIMOTU_R
      RcppThread::checkUserInterrupt();
#endif
    }
  }

} // namespace

HammingClusterWorker::HammingClusterWorker(
  const SequenceSet &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const int min_overlap,
  const bool ignore_gaps,
  int verbose,
  std::size_t worker_threads
) : DistClusterWorker(seq, clust_algo, pgb, verbose, worker_threads), pss(seq),
min_overlap(min_overlap), ignore_gaps(ignore_gaps) {
  tracked_pss_bytes = PackedSequenceSet::estimate_bytes(
    seq.size(), pss.alen
  );
  if (tracked_pss_bytes > 0) {
    if (auto *mca = dynamic_cast<MultipleClusterAlgorithm *>(&clust_algo)) {
      mca->acquire_memory(tracked_pss_bytes, "PackedSequenceSet");
    }
  }
}

HammingClusterWorker::~HammingClusterWorker() {
  if (tracked_pss_bytes > 0) {
    if (auto *mca = dynamic_cast<MultipleClusterAlgorithm *>(&clust_algo)) {
      mca->release_memory(tracked_pss_bytes);
    }
  }
}

template<int verbose>
HammingSplitClusterWorker<verbose>::HammingSplitClusterWorker(
  const SequenceSet &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const int min_overlap,
  const bool ignore_gaps,
  std::size_t worker_threads
) : HammingClusterWorker(seq, clust_algo, pgb, min_overlap, ignore_gaps,
   verbose, worker_threads) {};

template<int verbose>
void HammingSplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  std::size_t pairs_since_check = 0;

  for (size_t pg_index = begin; pg_index < pair_generators.size(); pg_index += threads) {
    auto & pg = pair_generators[pg_index];
    ClusterAlgorithm *my_algo = tile_algo(pg.get());
    if (!my_algo)
    {
      continue;
    }
    OPTIMOTU_DEBUG(
      2,
      << "HammingSplitClusterWorker thread " << pg_index
      << " entered" << std::endl
    );

    // Mapped tile children want tile-local indices; MCA / root want global.
    const bool use_local =
        dynamic_cast<MappedClusterAlgorithm *>(my_algo) != nullptr;

    auto consume = [&](std::size_t i0, std::size_t j0, j_t li, j_t lj)
    {
      ++my_prealigned;
      double d = pss.dist(
          static_cast<int>(i0),
          static_cast<int>(j0),
          min_overlap,
          ignore_gaps);
      if (d < 1.0)
        ++my_aligned;
      if (use_local)
      {
        double threshold = my_algo->max_relevant_local(li, lj);
        OPTIMOTU_DEBUG(
            4,
            << "thread" << pg_index
            << ": seqs " << lj << " (j0=" << j0 << ")"
            << " and " << li << " (i0=" << i0 << ")"
            << " max relevant=" << threshold
            << std::endl);
        OPTIMOTU_DEBUG(
            4,
            << (d <= threshold ? "*" : " ")
            << " distance=" << d
            << std::endl;);
        if (d <= threshold)
          my_algo->apply_local(li, lj, d);
      }
      else
      {
        double threshold = my_algo->max_relevant(
            static_cast<j_t>(i0), static_cast<j_t>(j0));
        OPTIMOTU_DEBUG(
            4,
            << "thread" << pg_index
            << ": seqs " << j0 << " and " << i0
            << " max relevant=" << threshold
            << std::endl);
        OPTIMOTU_DEBUG(
            4,
            << (d <= threshold ? "*" : " ")
            << " distance=" << d
            << std::endl;);
        if (d <= threshold)
          (*my_algo)(static_cast<j_t>(i0), static_cast<j_t>(j0), d);
      }
      hamming_maybe_interrupt(pairs_since_check);
    };

    if (auto *diag = dynamic_cast<AllPairGenerator *>(pg.get()))
    {
      const std::size_t offset = diag->seq_offset();
      const std::size_t n_tile = diag->seq_count();
      for (std::size_t i0 = offset + 1; i0 < offset + n_tile; ++i0)
      {
        const j_t li = static_cast<j_t>(i0 - offset);
        for (std::size_t j0 = offset; j0 < i0; ++j0)
        {
          consume(i0, j0, li, static_cast<j_t>(j0 - offset));
        }
      }
    }
    else if (auto *bip = dynamic_cast<BipartitePairGenerator *>(pg.get()))
    {
      const std::size_t begin_i = bip->seq_begin_i();
      const std::size_t end_i = bip->seq_end_i();
      const std::size_t begin_j = bip->seq_begin_j();
      const std::size_t end_j = bip->seq_end_j();
      const std::size_t nj = end_j - begin_j;
      for (std::size_t i0 = begin_i; i0 < end_i; ++i0)
      {
        const j_t li = static_cast<j_t>(i0 - begin_i + nj);
        for (std::size_t j0 = begin_j; j0 < end_j; ++j0)
        {
          consume(i0, j0, li, static_cast<j_t>(j0 - begin_j));
        }
      }
    }
    else
    {
      // Fallback for unexpected PairGenerator types.
      while (*pg)
      {
        std::size_t i0 = pg->i0();
        std::size_t j0 = pg->j0();
        consume(
            i0,
            j0,
            static_cast<j_t>(pg->i()),
            static_cast<j_t>(pg->j()));
        ++(*pg);
      }
    }

#ifdef OPTIMOTU_R
    RcppThread::checkUserInterrupt();
#endif
    pairs_since_check = 0;

    mutex.lock();
    OPTIMOTU_DEBUG(2, << "thread " << pg_index << " ready to merge" << std::endl);
    _aligned += my_aligned;
    _prealigned += my_prealigned;
    mutex.unlock();
    my_aligned = 0;
    my_prealigned = 0;
    finish_tile(my_algo);
    OPTIMOTU_DEBUG(2, << "thread " << pg_index << " done" << std::endl);
  }
}

template<int verbose>
HammingConcurrentClusterWorker<verbose>::HammingConcurrentClusterWorker(
  const SequenceSet &seq,
  ClusterAlgorithm &clust_algo,
  DivisiblePairGenerator::Builder & pgb,
  const int min_overlap,
  const bool ignore_gaps,
  std::size_t worker_threads
) : HammingClusterWorker(seq, clust_algo, pgb, min_overlap, ignore_gaps,
   verbose, worker_threads) {};

template<int verbose>
void HammingConcurrentClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  std::size_t pairs_since_check = 0;

  for (size_t pg_index = begin; pg_index < pair_generators.size(); pg_index += threads) {
    auto & pg = pair_generators[pg_index];
    OPTIMOTU_DEBUG(
      2,
      << "HammingConcurrentClusterWorker thread " << pg_index
      << " entered" << std::endl
    );

    auto consume = [&](std::size_t i0, std::size_t j0)
    {
      OPTIMOTU_DEBUG(
          4,
          << "thread" << pg_index
          << ": seqs " << j0 << " and " << i0
          << std::endl);
      double threshold = clust_algo.max_relevant(
          static_cast<j_t>(i0), static_cast<j_t>(j0));
      ++my_prealigned;
      double d = pss.dist(
          static_cast<int>(i0),
          static_cast<int>(j0),
          min_overlap,
          ignore_gaps);
      if (d < 1.0)
        ++my_aligned;
      if (d < threshold)
        clust_algo(static_cast<j_t>(i0), static_cast<j_t>(j0), d);
      hamming_maybe_interrupt(pairs_since_check);
    };

    if (auto *diag = dynamic_cast<AllPairGenerator *>(pg.get()))
    {
      const std::size_t offset = diag->seq_offset();
      const std::size_t n_tile = diag->seq_count();
      for (std::size_t i0 = offset + 1; i0 < offset + n_tile; ++i0)
      {
        for (std::size_t j0 = offset; j0 < i0; ++j0)
        {
          consume(i0, j0);
        }
      }
    }
    else if (auto *bip = dynamic_cast<BipartitePairGenerator *>(pg.get()))
    {
      const std::size_t begin_i = bip->seq_begin_i();
      const std::size_t end_i = bip->seq_end_i();
      const std::size_t begin_j = bip->seq_begin_j();
      const std::size_t end_j = bip->seq_end_j();
      for (std::size_t i0 = begin_i; i0 < end_i; ++i0)
      {
        for (std::size_t j0 = begin_j; j0 < end_j; ++j0)
        {
          consume(i0, j0);
        }
      }
    }
    else
    {
      while (*pg)
      {
        consume(pg->i0(), pg->j0());
        ++(*pg);
      }
    }

#ifdef OPTIMOTU_R
    RcppThread::checkUserInterrupt();
#endif
    pairs_since_check = 0;

    mutex.lock();
    _aligned += my_aligned;
    _prealigned += my_prealigned;
    OPTIMOTU_DEBUG(2, << "thread " << pg_index << " done" << std::endl);
    mutex.unlock();
    my_aligned = 0;
    my_prealigned = 0;
  }
}

template class HammingSplitClusterWorker<0>;
template class HammingSplitClusterWorker<1>;
template class HammingSplitClusterWorker<2>;
template class HammingSplitClusterWorker<3>;
template class HammingSplitClusterWorker<4>;
template class HammingConcurrentClusterWorker<0>;
template class HammingConcurrentClusterWorker<1>;
template class HammingConcurrentClusterWorker<2>;
template class HammingConcurrentClusterWorker<3>;
template class HammingConcurrentClusterWorker<4>;
