#include "SparseDistanceMatrix.h"
#include "optimotu.h"

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

void SparseDistanceMatrix::push_back(
    size_t s1, size_t s2,
    int sc1, int sc2,
    double d1, double d2, std::string cigar) {
  OPTIMOTU_STOP("attempted to add cigar to SparseDistanceMatrix");
}

void SparseDistanceMatrix::push_back(
    size_t s1, size_t s2,
    int sc1, int sc2,
    double d1, double d2, int al, int ni, int nd, int mi, int md) {
  OPTIMOTU_STOP("attempted to add gapstats to SparseDistanceMatrix");
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

void SparseDistanceMatrix::append(
    std::vector<size_t> &seq1new,
    std::vector<size_t> &seq2new,
    std::vector<int> &score1new,
    std::vector<int> &score2new,
    std::vector<double> &dist1new,
    std::vector<double> &dist2new,
    std::vector<int> &align_lengthnew,
    std::vector<int> &n_insertnew,
    std::vector<int> &n_deletenew,
    std::vector<int> &max_insertnew,
    std::vector<int> &max_deletenew) {
  OPTIMOTU_STOP("attempted to add gapstats to SparseDistanceMatrix");
}

void SparseDistanceMatrix::append(
    std::vector<size_t> &seq1new,
    std::vector<size_t> &seq2new,
    std::vector<int> &score1new,
    std::vector<int> &score2new,
    std::vector<double> &dist1new,
    std::vector<double> &dist2new,
    std::vector<std::string> &cigarnew) {
  OPTIMOTU_STOP("attempted to add cigar to SparseDistanceMatrix");
}

SparseDistanceMatrixCigar::SparseDistanceMatrixCigar(
    std::vector<size_t> &seq1,
    std::vector<size_t> &seq2,
    std::vector<int> &score1,
    std::vector<int> &score2,
    std::vector<double> &dist1,
    std::vector<double> &dist2,
    std::vector<std::string> &cigar
) : SparseDistanceMatrix(seq1, seq2, score1, score2, dist1, dist2),
   cigar(cigar)
   {}

void SparseDistanceMatrixCigar::push_back(
  size_t s1, size_t s2,
  int sc1, int sc2,
  double d1, double d2) {
    OPTIMOTU_STOP("attempted to push to SparseDistanceMatrixCigar without cigar");
}

void SparseDistanceMatrixCigar::push_back(
    size_t s1, size_t s2,
    int sc1, int sc2,
    double d1, double d2, std::string cigar) {
  SparseDistanceMatrix::push_back(s1, s2, sc1, sc2, d1, d2);
  this->cigar.push_back(cigar);
}

void SparseDistanceMatrixCigar::append(
    std::vector<size_t> &seq1new,
    std::vector<size_t> &seq2new,
    std::vector<int> &score1new,
    std::vector<int> &score2new,
    std::vector<double> &dist1new,
    std::vector<double> &dist2new) {
  OPTIMOTU_STOP("attempted to append to SparseDistanceMatrixCigar without cigar");
}

void SparseDistanceMatrixCigar::append(
    std::vector<size_t> &seq1new,
    std::vector<size_t> &seq2new,
    std::vector<int> &score1new,
    std::vector<int> &score2new,
    std::vector<double> &dist1new,
    std::vector<double> &dist2new,
    std::vector<std::string> &cigarnew) {
  SparseDistanceMatrix::append(seq1new, seq2new, score1new, score2new, dist1new, dist2new);
  this->cigar.insert(this->cigar.end(), cigarnew.begin(), cigarnew.end());
  cigarnew.clear();
}

SparseDistanceMatrixGapstats::SparseDistanceMatrixGapstats(
    std::vector<size_t> &seq1,
    std::vector<size_t> &seq2,
    std::vector<int> &score1,
    std::vector<int> &score2,
    std::vector<double> &dist1,
    std::vector<double> &dist2,
    std::vector<int> &align_length,
    std::vector<int> &n_insert,
    std::vector<int> &n_delete,
    std::vector<int> &max_insert,
    std::vector<int> &max_delete
) : SparseDistanceMatrix(seq1, seq2, score1, score2, dist1, dist2),
  align_length(align_length),
  n_insert(n_insert),
  n_delete(n_delete),
  max_insert(max_insert),
  max_delete(max_delete) {}

void SparseDistanceMatrixGapstats::push_back(
    size_t s1, size_t s2,
    int sc1, int sc2,
    double d1, double d2) {
  OPTIMOTU_STOP("attempted to push to SparseDistanceMatrixGapstats without gapstats");
}

void SparseDistanceMatrixGapstats::push_back(
    size_t s1, size_t s2,
    int sc1, int sc2,
    double d1, double d2, int al, int ni, int nd, int mi, int md) {
  SparseDistanceMatrix::push_back(s1, s2, sc1, sc2, d1, d2);
  this->align_length.push_back(al);
  this->n_insert.push_back(ni);
  this->n_delete.push_back(nd);
  this->max_insert.push_back(mi);
  this->max_delete.push_back(md);
}

void SparseDistanceMatrixGapstats::append(
    std::vector<size_t> &seq1new,
    std::vector<size_t> &seq2new,
    std::vector<int> &score1new,
    std::vector<int> &score2new,
    std::vector<double> &dist1new,
    std::vector<double> &dist2new) {
  OPTIMOTU_STOP("attempted to append to SparseDistanceMatrixGapstats without gapstats");
}

void SparseDistanceMatrixGapstats::append(
    std::vector<size_t> &seq1new,
    std::vector<size_t> &seq2new,
    std::vector<int> &score1new,
    std::vector<int> &score2new,
    std::vector<double> &dist1new,
    std::vector<double> &dist2new,
    std::vector<int> &align_lengthnew,
    std::vector<int> &n_insertnew,
    std::vector<int> &n_deletenew,
    std::vector<int> &max_insertnew,
    std::vector<int> &max_deletenew) {
  SparseDistanceMatrix::append(seq1new, seq2new, score1new, score2new, dist1new, dist2new);
  this->align_length.insert(this->align_length.end(), align_lengthnew.begin(), align_lengthnew.end());
  align_lengthnew.clear();
  this->n_insert.insert(this->n_insert.end(), n_insertnew.begin(), n_insertnew.end());
  n_insertnew.clear();
  this->n_delete.insert(this->n_delete.end(), n_deletenew.begin(), n_deletenew.end());
  n_deletenew.clear();
  this->max_insert.insert(this->max_insert.end(), max_insertnew.begin(), max_insertnew.end());
  max_insertnew.clear();
  this->max_delete.insert(this->max_delete.end(), max_deletenew.begin(), max_deletenew.end());
  max_deletenew.clear();
}