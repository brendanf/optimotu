#ifndef OPTIMOTU_SPARSEDISTANCEMATRIX_H_INCLUDED
#define OPTIMOTU_SPARSEDISTANCEMATRIX_H_INCLUDED

#include <vector>
#include <mutex>
#include <string>

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

  virtual void push_back(size_t s1, size_t s2, int sc1, int sc2, double d1, double d2);

  virtual void push_back(size_t s1, size_t s2, int sc1, int sc2, double d1, double d2, std::string cigar);

  virtual void push_back(size_t s1, size_t s2, int sc1, int sc2, double d1, double d2, int al, int ni, int nd, int mi, int md);

  virtual void append(
      std::vector<size_t> &seq1new,
      std::vector<size_t> &seq2new,
      std::vector<int> &score1new,
      std::vector<int> &score2new,
      std::vector<double> &dist1new,
      std::vector<double> &dist2new
  );

  virtual void append(
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
    std::vector<int> &max_deletenew
  );

  virtual void append(
    std::vector<size_t> &seq1new,
    std::vector<size_t> &seq2new,
    std::vector<int> &score1new,
    std::vector<int> &score2new,
    std::vector<double> &dist1new,
    std::vector<double> &dist2new,
    std::vector<std::string> &cigarnew
  );

  virtual ~SparseDistanceMatrix() = default;
};

struct SparseDistanceMatrixCigar : public SparseDistanceMatrix {
  std::vector<std::string> &cigar;

  SparseDistanceMatrixCigar(
    std::vector<size_t> &seq1,
    std::vector<size_t> &seq2,
    std::vector<int> &score1,
    std::vector<int> &score2,
    std::vector<double> &dist1,
    std::vector<double> &dist2,
    std::vector<std::string> &cigar
  );

  virtual void push_back(size_t s1, size_t s2, int sc1, int sc2, double d1, double d2) override;

  virtual void push_back(size_t s1, size_t s2, int sc1, int sc2, double d1, double d2, std::string cigar) override;

  virtual void append(
    std::vector<size_t> &seq1new,
    std::vector<size_t> &seq2new,
    std::vector<int> &score1new,
    std::vector<int> &score2new,
    std::vector<double> &dist1new,
    std::vector<double> &dist2new
  ) override;

  virtual void append(
    std::vector<size_t> &seq1new,
    std::vector<size_t> &seq2new,
    std::vector<int> &score1new,
    std::vector<int> &score2new,
    std::vector<double> &dist1new,
    std::vector<double> &dist2new,
    std::vector<std::string> &cigarnew
  ) override;

  virtual ~SparseDistanceMatrixCigar() = default;
};

struct SparseDistanceMatrixGapstats : public SparseDistanceMatrix {
  std::vector<int> &align_length, &n_insert, &n_delete, &max_insert, &max_delete;

  SparseDistanceMatrixGapstats(
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
  );

  virtual void push_back(size_t s1, size_t s2, int sc1, int sc2, double d1, double d2) override;

  virtual void push_back(size_t s1, size_t s2, int sc1, int sc2, double d1, double d2, int al, int ni, int nd, int mi, int md) override;

  virtual void append(
    std::vector<size_t> &seq1new,
    std::vector<size_t> &seq2new,
    std::vector<int> &score1new,
    std::vector<int> &score2new,
    std::vector<double> &dist1new,
    std::vector<double> &dist2new
  ) override;

  virtual void append(
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
    std::vector<int> &max_deletenew
  ) override;

  virtual ~SparseDistanceMatrixGapstats() = default;
        };

#endif //OPTIMOTU_SPARSEDISTANCEMATRIX_H_INCLUDED
