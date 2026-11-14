// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXMODPRIME_H_
#define SRC_MATRIX_MAT_MATRIXMODPRIME_H_

// clang-format off
#include "ResidueArithmetic.h"
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

/*
  Linear algebra modulo a prime p <= 2^62 on matrices of residues stored
  as rows of int64_t, through the residue operations of
  ResidueArithmetic.h (no product overflows). The entries given may be any
  integers; they are reduced first.
 */

using ModMatrix = std::vector<std::vector<int64_t>>;

// X Y modulo mod.
inline ModMatrix ModMatrixProduct(ModMatrix const &X, ModMatrix const &Y,
                                  int64_t mod) {
  size_t n = X.size();
  size_t k = Y.size();
  size_t l = Y.empty() ? 0 : Y[0].size();
  ModMatrix Z(n, std::vector<int64_t>(l, 0));
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < k; j++) {
      int64_t x = ResMod(X[i][j], mod);
      if (x == 0) {
        continue;
      }
      for (size_t u = 0; u < l; u++) {
        Z[i][u] = ResAdd(Z[i][u], ResMul(x, ResMod(Y[j][u], mod), mod), mod);
      }
    }
  }
  return Z;
}

inline ModMatrix ModMatrixIdentity(size_t n) {
  ModMatrix Id(n, std::vector<int64_t>(n, 0));
  for (size_t i = 0; i < n; i++) {
    Id[i][i] = 1;
  }
  return Id;
}

/*
  The reduced row echelon form of the row space modulo the prime p: the
  basis rows and their pivot columns. The coordinates of a vector of the
  row space in that basis are its entries at the pivots.
 */
struct RowSpaceModP {
  ModMatrix basis;
  std::vector<size_t> pivots;
};

inline RowSpaceModP ComputeRowSpaceModP(ModMatrix rows, int64_t p) {
  RowSpaceModP res;
  size_t n_col = rows.empty() ? 0 : rows[0].size();
  size_t r = 0;
  for (auto &row : rows) {
    for (auto &x : row) {
      x = ResMod(x, p);
    }
  }
  for (size_t col = 0; col < n_col && r < rows.size(); col++) {
    size_t piv = r;
    while (piv < rows.size() && rows[piv][col] == 0) {
      piv++;
    }
    if (piv == rows.size()) {
      continue;
    }
    std::swap(rows[r], rows[piv]);
    int64_t inv = ResInv(rows[r][col], p);
    for (size_t j = col; j < n_col; j++) {
      rows[r][j] = ResMul(rows[r][j], inv, p);
    }
    for (size_t i = 0; i < rows.size(); i++) {
      int64_t c = rows[i][col];
      if (i != r && c != 0) {
        for (size_t j = col; j < n_col; j++) {
          if (rows[r][j] != 0) {
            rows[i][j] = ResSub(rows[i][j], ResMul(c, rows[r][j], p), p);
          }
        }
      }
    }
    res.pivots.push_back(col);
    r++;
  }
  rows.resize(r);
  res.basis = std::move(rows);
  return res;
}

// The rank modulo the prime p.
inline size_t RankModPrime(ModMatrix const &rows, int64_t p) {
  return ComputeRowSpaceModP(rows, p).pivots.size();
}

// The coordinates of a vector of the row space in its echelon basis.
inline std::vector<int64_t>
CoordinatesInRowSpaceModP(RowSpaceModP const &rs,
                          std::vector<int64_t> const &v) {
  std::vector<int64_t> c(rs.pivots.size());
  for (size_t i = 0; i < rs.pivots.size(); i++) {
    c[i] = v[rs.pivots[i]];
  }
  return c;
}

// A basis of the x with x M = 0 modulo p (M with r rows).
inline ModMatrix LeftNullspaceModP(ModMatrix const &M, int64_t p) {
  size_t r = M.size();
  size_t c = M.empty() ? 0 : M[0].size();
  ModMatrix aug(r, std::vector<int64_t>(c + r, 0));
  for (size_t i = 0; i < r; i++) {
    for (size_t j = 0; j < c; j++) {
      aug[i][j] = ResMod(M[i][j], p);
    }
    aug[i][c + i] = 1;
  }
  RowSpaceModP rs = ComputeRowSpaceModP(std::move(aug), p);
  ModMatrix res;
  for (size_t i = 0; i < rs.basis.size(); i++) {
    if (rs.pivots[i] >= c) {
      res.push_back(std::vector<int64_t>(rs.basis[i].begin() + c,
                                         rs.basis[i].end()));
    }
  }
  return res;
}

// The solution x of x B = v modulo p (B with independent rows), if any.
inline std::optional<std::vector<int64_t>>
SolveRowsModP(ModMatrix const &B, std::vector<int64_t> const &v, int64_t p) {
  size_t r = B.size();
  ModMatrix M = B;
  M.push_back(v);
  ModMatrix ker = LeftNullspaceModP(M, p);
  for (auto const &x : ker) {
    if (x[r] != 0) {
      int64_t inv = ResInv(ResNeg(x[r], p), p);
      std::vector<int64_t> sol(r);
      for (size_t i = 0; i < r; i++) {
        sol[i] = ResMul(x[i], inv, p);
      }
      return sol;
    }
  }
  return {};
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXMODPRIME_H_
// clang-format on
