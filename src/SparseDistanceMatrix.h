#ifndef _SparseDistanceMatrix_
#define _SparseDistanceMatrix_
#include <vector>
#include <mutex>

struct SparseDistanceMatrix {
  std::vector<size_t> &seq1, &seq2;
  std::vector<int> &score1, &score2;
  std::vector<double> &dist1, &dist2;
  std::mutex mutex;

  SparseDistanceMatrix(
    std::vector<size_t> &seq1,
    std::vector<size_t> &seq2,
    std::vector<int> &score1,
    std::vector<int> &score2,
    std::vector<double> &dist1,
    std::vector<double> &dist2
  );

  void push_back(size_t s1, size_t s2, int sc1, int sc2, double d1, double d2);

  void append(
      std::vector<size_t> &seq1new,
      std::vector<size_t> &seq2new,
      std::vector<int> &score1new,
      std::vector<int> &score2new,
      std::vector<double> &dist1new,
      std::vector<double> &dist2new
  );
};
#endif
