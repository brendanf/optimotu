#include "HammingClusterWorker.h"

typedef RcppParallel::RMatrix<int> matrix_t;

PackedSequenceSet::PackedSequenceSet(const std::vector<std::string> &seq) {
  num_seqs = seq.size();
  if (num_seqs == 0) return;
  packed_seq.reserve(seq.size());
  mask.reserve(seq.size());

  alen = seq[0].size();
  ulen = alen / NUCLEOTIDES_IN_WORD;
  if (alen > ulen * NUCLEOTIDES_IN_WORD)
    ulen++;
  mulen = alen / NUCLEOTIDES_IN_WORD / 4;
  if (alen > mulen * NUCLEOTIDES_IN_WORD / 4)
    mulen++;
  for (auto s : seq) {
    if (s.size() != alen)
      OPTIMOTU_STOP("PackedSequenceSet: All sequences must have the same length.\n");
    packed_seq.emplace_back(ulen);
    mask.emplace_back(mulen);
    nucleotide2binary(s.c_str(), alen, &packed_seq.back()[0], &mask.back()[0]);
  }
}

double PackedSequenceSet::dist(const int i, const int j) const {
  return pdistB(&packed_seq[i][0], &mask[i][0], &packed_seq[j][0], &mask[j][0], ulen, mulen);
}

HammingClusterWorker::HammingClusterWorker(
  const std::vector<std::string> &seq,
  ClusterAlgorithm &clust_algo,
  const uint8_t threads,
  bool verbose
) : pss(seq), clust_algo(clust_algo),
threads(threads), verbose(verbose) {};

HammingSplitClusterWorker::HammingSplitClusterWorker(
  const std::vector<std::string> &seq,
  ClusterAlgorithm &clust_algo,
  const uint8_t threads,
  bool verbose
) : HammingClusterWorker(seq, clust_algo, threads, false) {};

void HammingSplitClusterWorker::operator()(std::size_t begin, std::size_t end) {
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
  if (verbose) {
    mutex.lock();
    std::cout << "HammingSplitClusterWorker thread " << begin << " entered; sequences [" <<
      begin_i << ", "<< end_i << ")" << std::endl;
    mutex.unlock();
  }
  for (size_t i = begin_i; i < end_i; i++) {
    for (size_t j = 0; j < i; j++) {
      double threshold = my_algo->max_relevant(i, j);
      ++my_prealigned;
      double d = pss.dist(i, j);
      if (d < 1.0) ++my_aligned;
      if (d < threshold) (*my_algo)(j, i, d);
      RcppThread::checkUserInterrupt();
    }
  }
  mutex.lock();
  if (verbose) std::cout << "thread " << begin << " ready to merge" << std::endl;
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  mutex.unlock();
  my_algo->merge_into_parent();
  if (verbose) std::cout << "thread " << begin << " done" << std::endl;
}

void HammingConcurrentClusterWorker::operator()(std::size_t begin, std::size_t end) {
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
  if (verbose) {
    mutex.lock();
    std::cout << "HammingConcurrentClusterWorker thread " << begin
              << " entered; sequences [" << begin_i
              << ", "<< end_i << ")" << std::endl;
    mutex.unlock();
  }
  for (size_t i = begin_i; i < end_i; i++) {
    for (size_t j = 0; j < i; j++) {
      // mutex.lock();
      // std::cout << "Thread " << begin
      //           << ": seqs " << j
      //           << " and " << i
      //           << std::endl;
      // mutex.unlock();
      double threshold = clust_algo.max_relevant(i, j);
      ++my_prealigned;
      double d = pss.dist(i, j);
      if (d < 1.0) ++my_aligned;
      if (d < threshold) clust_algo(j, i, d);
      RcppThread::checkUserInterrupt();
    }
  }
  mutex.lock();
  _aligned += my_aligned;
  _prealigned += my_prealigned;
  // std::cout << "Exiting thread " << begin << std::endl;
  mutex.unlock();
}
