// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_Matrix.h"
#include "MAT_SparseMatrixMod.h"
#include <string>
#include <vector>
// clang-format on

/*
  Tests of the sparse linear algebra modulo p against the dense rational
  rank (equal to the rank modulo p for random matrices and a large p, up
  to a probability of order 1/p):
  - the number of independent rows found by SparseEchelonMod is the rank,
  - SparseLeftKernelMod returns n_row - rank vectors, all in the kernel,
  - the normal form of v + (combination of the rows) is the one of v.
 */

int64_t const p_test = 2147483647;

SparseVectorMod random_sparse_row(int n_col, int n_nonzero) {
  std::map<int, int64_t> map;
  for (int k = 0; k < n_nonzero; k++) {
    int idx = random_int(0, n_col - 1);
    map[idx] = ResidueModPrime(random_int(-3, 3), p_test);
  }
  SparseVectorMod v;
  for (auto &kv : map) {
    if (kv.second != 0) {
      v.push_back(kv);
    }
  }
  return v;
}

MyMatrix<mpq_class> to_dense(std::vector<SparseVectorMod> const &l_row, int n_col) {
  int n_row = l_row.size();
  MyMatrix<mpq_class> M = ZeroMatrix<mpq_class>(n_row, n_col);
  for (int i = 0; i < n_row; i++) {
    for (auto &e : l_row[i]) {
      int64_t val = e.second;
      if (val > p_test / 2) {
        val -= p_test;
      }
      M(i, e.first) = mpq_class(static_cast<long>(val));
    }
  }
  return M;
}

bool test_one(int n_row, int n_col, int n_nonzero) {
  std::vector<SparseVectorMod> l_row;
  for (int i = 0; i < n_row; i++) {
    l_row.push_back(random_sparse_row(n_col, n_nonzero));
  }
  // Some dependent rows.
  for (int k = 0; k < n_row / 4; k++) {
    int i1 = random_int(0, n_row - 1);
    int i2 = random_int(0, n_row - 1);
    l_row.push_back(SparseVectorAddMultipleMod(l_row[i1], random_int(1, 5), l_row[i2], p_test));
  }
  int n_row_tot = l_row.size();
  int rank_dense = RankMat(to_dense(l_row, n_col));
  SparseEchelonMod echelon(p_test);
  for (auto &v : l_row) {
    echelon.insert(v);
  }
  if (static_cast<int>(echelon.rank()) != rank_dense) {
    std::cerr << "The sparse rank " << echelon.rank() << " differs from the dense rank " << rank_dense << "\n";
    return false;
  }
  std::vector<SparseVectorMod> l_kernel = SparseLeftKernelMod(l_row, n_col, p_test);
  if (static_cast<int>(l_kernel.size()) != n_row_tot - rank_dense) {
    std::cerr << "The kernel has dimension " << l_kernel.size() << " instead of " << n_row_tot - rank_dense << "\n";
    return false;
  }
  for (auto &x : l_kernel) {
    SparseVectorMod sum;
    for (auto &e : x) {
      sum = SparseVectorAddMultipleMod(sum, e.second, l_row[e.first], p_test);
    }
    if (!sum.empty()) {
      std::cerr << "A kernel vector is not in the kernel\n";
      return false;
    }
  }
  for (int k = 0; k < 10; k++) {
    SparseVectorMod v = random_sparse_row(n_col, n_nonzero);
    SparseVectorMod w = v;
    for (int u = 0; u < 3; u++) {
      int i = random_int(0, n_row_tot - 1);
      w = SparseVectorAddMultipleMod(w, random_int(-7, 7), l_row[i], p_test);
    }
    if (echelon.reduce(v) != echelon.reduce(w)) {
      std::cerr << "The normal form depends on the representative\n";
      return false;
    }
  }
  return true;
}

int main() {
  HumanTime time;
  try {
    int n_error = 0;
    for (int iter = 0; iter < 50; iter++) {
      int n_row = random_int(1, 60);
      int n_col = random_int(1, 60);
      int n_nonzero = random_int(1, 6);
      if (!test_one(n_row, n_col, n_nonzero)) {
        std::cerr << "Failure at iter=" << iter << " n_row=" << n_row << " n_col=" << n_col << "\n";
        n_error += 1;
      }
    }
    if (n_error > 0) {
      std::cerr << "Test_SparseMatrixMod: " << n_error << " failures\n";
      return 1;
    }
    std::cerr << "Test_SparseMatrixMod: all tests passed\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_SparseMatrixMod\n";
    exit(e.eVal);
  }
  runtime(time);
}
