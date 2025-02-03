#include "HammingClusterWorker.h"
extern "C" {
#include "defs.h"
}

typedef RcppParallel::RMatrix<int> matrix_t;

PackedSequenceSet::PackedSequenceSet(const std::vector<std::string> &seq) {
  num_seqs = seq.size();
  if (num_seqs == 0) return;
  packed_seq.reserve(seq.size());
  mask.reserve(seq.size());
  start.reserve(seq.size());
  end.reserve(seq.size());

  for (auto s : seq) {
    // if (s.size() != alen)
    //   OPTIMOTU_STOP("PackedSequenceSet: All sequences must have the same length.\n");
    alen = s.size();
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
    nucleotide2binary(s.c_str(), alen, &packed_seq.back()[0], &mask.back()[0], &start.back(), &end.back());
  }
}

double PackedSequenceSet::dist(const int i, const int j, const int min_overlap, const bool ignore_gap) const {
  int s, e;
  s = start[j] >= start[i] ? start[j] : start[i];
  e = end[j] <= end[i] ? end[j] : end[i];
  if (s >= e) return 1.0;
  return pdistB(&packed_seq[i][0], &mask[i][0],
                &packed_seq[j][0], &mask[j][0],
                s, e, min_overlap);
}


HammingClusterWorker::HammingClusterWorker(
  const std::vector<std::string> &seq,
  ClusterAlgorithm &clust_algo,
  const std::uint8_t threads,
  const int min_overlap,
  const bool ignore_gaps
) : AlignClusterWorker(seq, clust_algo, threads), pss(seq),
min_overlap(min_overlap), ignore_gaps(ignore_gaps) {};

template<int verbose>
HammingSplitClusterWorker<verbose>::HammingSplitClusterWorker(
  const std::vector<std::string> &seq,
  ClusterAlgorithm &clust_algo,
  const std::uint8_t threads,
  const int min_overlap,
  const bool ignore_gaps
) : HammingClusterWorker(seq, clust_algo, threads, min_overlap, ignore_gaps) {};

template<int verbose>
void HammingSplitClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {
  double n = pss.num_seqs;
  double m = (n*n - 3.0*n + 2.0)/2.0;
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  size_t begin_i;

  ClusterAlgorithm * my_algo = clust_algo.make_child();

  if (begin == 0) {
    begin_i = 1;
  } else {
    begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
  }
  size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
  OPTIMOTU_DEBUG(
    1,
    << "HammingSplitClusterWorker thread " << begin
    << " entered; sequences [" << begin_i
    << ", "<< end_i << ")" << std::endl
  );
  for (size_t i = begin_i; i < end_i; i++) {
    for (size_t j = 0; j < i; j++) {
      double threshold = my_algo->max_relevant(i, j);
      OPTIMOTU_DEBUG(
        2,
        << "thread" << begin
        << ": seqs " << j
        << " and " << i
        << " max relevant=" << threshold
        << std::endl
      );
      ++my_prealigned;
      double d = pss.dist(i, j, min_overlap, ignore_gaps);
      if (d < 1.0) ++my_aligned;

      OPTIMOTU_DEBUG(
        2,
        << (d <= threshold ? "*" : " ")
        << " distance=" << d
        << std::endl;
      );
      if (d <= threshold) (*my_algo)(j, i, d);
      RcppThread::checkUserInterrupt();
    }
  }
  mutex.lock();
  OPTIMOTU_DEBUG(1, << "thread " << begin << " ready to merge" << std::endl);
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  mutex.unlock();
  my_algo->merge_into_parent();
  OPTIMOTU_DEBUG(1, << "thread " << begin << " done" << std::endl);
}

template<int verbose>
HammingConcurrentClusterWorker<verbose>::HammingConcurrentClusterWorker(
  const std::vector<std::string> &seq,
  ClusterAlgorithm &clust_algo,
  const std::uint8_t threads,
  const int min_overlap,
  const bool ignore_gaps
) : HammingClusterWorker(seq, clust_algo, threads, min_overlap, ignore_gaps) {};

template<int verbose>
void HammingConcurrentClusterWorker<verbose>::operator()(std::size_t begin, std::size_t end) {

  double n = pss.num_seqs;
  double m = (n*n - 3.0*n + 2.0)/2.0;
  size_t my_prealigned = 0;
  size_t my_aligned = 0;
  size_t begin_i;

  if (begin == 0) {
    begin_i = 1;
  } else {
    begin_i = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*begin)/threads - 1.0)));
  }
  size_t end_i   = round(1.5 + 0.5*sqrt(9.0 + 8.0*((m*end)/threads - 1.0)));
  OPTIMOTU_DEBUG(
    1,
    << "HammingConcurrentClusterWorker thread " << begin
    << " entered; sequences [" << begin_i
    << ", "<< end_i << ")" << std::endl
  );
  for (size_t i = begin_i; i < end_i; i++) {
    for (size_t j = 0; j < i; j++) {
      OPTIMOTU_DEBUG(
        2,
        << "thread" << begin
        << ": seqs " << j
        << " and " << i
        << std::endl
      );

      double threshold = clust_algo.max_relevant(i, j);
      ++my_prealigned;
      double d = pss.dist(i, j, min_overlap, ignore_gaps);
      if (d < 1.0) ++my_aligned;
      if (d < threshold) clust_algo(j, i, d);
      RcppThread::checkUserInterrupt();
    }
  }
  mutex.lock();
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  OPTIMOTU_DEBUG(1, << "Exiting thread " << begin << std::endl);
  mutex.unlock();
}

template class HammingSplitClusterWorker<0>;
template class HammingSplitClusterWorker<1>;
template class HammingSplitClusterWorker<2>;
template class HammingConcurrentClusterWorker<0>;
template class HammingConcurrentClusterWorker<1>;
template class HammingConcurrentClusterWorker<2>;
