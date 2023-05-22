#include "SparseDistanceMatrix.h"

SparseDistanceMatrix::SparseDistanceMatrix(
  std::vector<size_t> &seq1,
  std::vector<size_t> &seq2,
  std::vector<int> &score1,
  std::vector<int> &score2,
  std::vector<double> &dist1,
  std::vector<double> &dist2
) : seq1(seq1), seq2(seq2),
score1(score1), score2(score2),
dist1(dist1), dist2(dist2),
mutex() {};

void SparseDistanceMatrix::push_back(
    size_t s1, size_t s2,
    int sc1, int sc2,
    double d1, double d2) {
  seq1.push_back(s1);
  seq2.push_back(s2);
  score1.push_back(sc1);
  score2.push_back(sc2);
  dist1.push_back(d1);
  dist2.push_back(d2);
}

void SparseDistanceMatrix::append(
    std::vector<size_t> &seq1new,
    std::vector<size_t> &seq2new,
    std::vector<int> &score1new,
    std::vector<int> &score2new,
    std::vector<double> &dist1new,
    std::vector<double> &dist2new
) {
  mutex.lock();
  std::copy(seq1new.begin(), seq1new.end(), std::back_inserter(seq1));
  seq1new.clear();
  std::copy(seq2new.begin(), seq2new.end(), std::back_inserter(seq2));
  seq2new.clear();
  std::copy(score1new.begin(), score1new.end(), std::back_inserter(score1));
  score1new.clear();
  std::copy(score2new.begin(), score2new.end(), std::back_inserter(score2));
  score2new.clear();
  std::copy(dist1new.begin(), dist1new.end(), std::back_inserter(dist1));
  dist1new.clear();
  std::copy(dist2new.begin(), dist2new.end(), std::back_inserter(dist2));
  dist2new.clear();
  mutex.unlock();
}
