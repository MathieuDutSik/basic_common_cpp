// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXMODPRIME_H_
#define SRC_MATRIX_MAT_MATRIXMODPRIME_H_

// clang-format off
#include "MAT_Matrix.h"
#include "ResidueArithmetic.h"
#include <cstdint>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

/*
  Linear algebra modulo a prime p <= 2^62 (and the products modulo any
  modulus) on matrices of residues in [0, m), through the residue
  operations of ResidueArithmetic.h, so that no product overflows.
 */

using ModMatrix = MyMatrix<T_residue>;
using ModVector = std::vector<T_residue>;

/*
  The residues modulo mod of a matrix with integer entries: T an integer
  type, or a field type whose entries have to be integers.
 */
template <typename T>
ModMatrix ModMatrixOfMatrix(MyMatrix<T> const &M, T_residue mod) {
  int n = M.rows();
  int l = M.cols();
  ModMatrix R(n, l);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < l; j++) {
      if constexpr (is_ring_field<T>::value) {
        if (!IsInteger(M(i, j))) {
          std::cerr << "MAT_MOD: ModMatrixOfMatrix of a non-integral entry\n";
          throw TerminalException{1};
        }
        R(i, j) = ResModOf(GetNumerator_z(M(i, j)), mod);
      } else {
        R(i, j) = ResModOf(M(i, j), mod);
      }
    }
  }
  return R;
}

// X Y modulo mod.
inline ModMatrix ModMatrixProduct(ModMatrix const &X, ModMatrix const &Y,
                                  T_residue mod) {
  int n = X.rows();
  int k = Y.rows();
  int l = Y.cols();
  ModMatrix Z = ZeroMatrix<T_residue>(n, l);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < k; j++) {
      T_residue x = X(i, j) % mod;
      if (x != 0) {
        for (int u = 0; u < l; u++) {
          Z(i, u) = ResAdd(Z(i, u), ResMul(x, Y(j, u) % mod, mod), mod);
        }
      }
    }
  }
  return Z;
}

inline ModMatrix ModMatrixIdentity(int n) { return IdentityMat<T_residue>(n); }

// The entries of X row after row.
inline ModVector FlattenModMatrix(ModMatrix const &X) {
  ModVector v;
  v.reserve(static_cast<size_t>(X.rows()) * static_cast<size_t>(X.cols()));
  for (int i = 0; i < X.rows(); i++) {
    for (int j = 0; j < X.cols(); j++) {
      v.push_back(X(i, j));
    }
  }
  return v;
}

// X with the row v appended (X may be empty, of any number of columns).
inline void AppendModRow(ModMatrix &X, ModVector const &v) {
  int r = X.rows();
  X.conservativeResize(r + 1, static_cast<int>(v.size()));
  for (size_t j = 0; j < v.size(); j++) {
    X(r, static_cast<int>(j)) = v[j];
  }
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

inline RowSpaceModP ComputeRowSpaceModP(ModMatrix rows, T_residue p) {
  RowSpaceModP res;
  int n_row = rows.rows();
  int n_col = rows.cols();
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++) {
      rows(i, j) %= p;
    }
  }
  int r = 0;
  for (int col = 0; col < n_col && r < n_row; col++) {
    int piv = r;
    while (piv < n_row && rows(piv, col) == 0) {
      piv++;
    }
    if (piv == n_row) {
      continue;
    }
    if (piv != r) {
      rows.row(r).swap(rows.row(piv));
    }
    T_residue inv = ResInv(rows(r, col), p);
    for (int j = col; j < n_col; j++) {
      rows(r, j) = ResMul(rows(r, j), inv, p);
    }
    for (int i = 0; i < n_row; i++) {
      T_residue c = rows(i, col);
      if (i != r && c != 0) {
        for (int j = col; j < n_col; j++) {
          if (rows(r, j) != 0) {
            rows(i, j) = ResSub(rows(i, j), ResMul(c, rows(r, j), p), p);
          }
        }
      }
    }
    res.pivots.push_back(static_cast<size_t>(col));
    r++;
  }
  res.basis = rows.topRows(r);
  return res;
}

/*
  The kernel of the system of the rows of rs (n_var unknowns), on the
  basis of the free unknowns: the vector of the free unknown f has a 1 at
  f, 0 at the other free unknowns, and minus the column f of the echelon
  form at the pivots.
 */
inline ModMatrix FreeKernelBasisModP(RowSpaceModP const &rs, int n_var,
                                     T_residue p) {
  std::vector<bool> is_pivot(n_var, false);
  for (size_t col : rs.pivots) {
    if (col < static_cast<size_t>(n_var)) {
      is_pivot[col] = true;
    }
  }
  std::vector<int> free_col;
  for (int col = 0; col < n_var; col++) {
    if (!is_pivot[col]) {
      free_col.push_back(col);
    }
  }
  int t = static_cast<int>(free_col.size());
  ModMatrix basis = ZeroMatrix<T_residue>(t, n_var);
  for (int f = 0; f < t; f++) {
    basis(f, free_col[f]) = 1;
    for (size_t r = 0; r < rs.pivots.size(); r++) {
      if (rs.pivots[r] < static_cast<size_t>(n_var)) {
        basis(f, static_cast<int>(rs.pivots[r])) =
            ResNeg(rs.basis(static_cast<int>(r), free_col[f]), p);
      }
    }
  }
  return basis;
}

// The rank modulo the prime p.
inline size_t RankModPrime(ModMatrix const &rows, T_residue p) {
  return ComputeRowSpaceModP(rows, p).pivots.size();
}

// The coordinates of a vector of the row space in its echelon basis.
inline ModVector CoordinatesInRowSpaceModP(RowSpaceModP const &rs,
                                           ModVector const &v) {
  ModVector c(rs.pivots.size());
  for (size_t i = 0; i < rs.pivots.size(); i++) {
    c[i] = v[rs.pivots[i]];
  }
  return c;
}

// A basis of the x with x M = 0 modulo p (M with r rows).
inline ModMatrix LeftNullspaceModP(ModMatrix const &M, T_residue p) {
  int r = M.rows();
  int c = M.cols();
  ModMatrix aug = ZeroMatrix<T_residue>(r, c + r);
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < c; j++) {
      aug(i, j) = M(i, j) % p;
    }
    aug(i, c + i) = 1;
  }
  RowSpaceModP rs = ComputeRowSpaceModP(std::move(aug), p);
  std::vector<int> kept;
  for (size_t i = 0; i < rs.pivots.size(); i++) {
    if (rs.pivots[i] >= static_cast<size_t>(c)) {
      kept.push_back(static_cast<int>(i));
    }
  }
  ModMatrix res(static_cast<int>(kept.size()), r);
  for (size_t u = 0; u < kept.size(); u++) {
    res.row(static_cast<int>(u)) = rs.basis.row(kept[u]).tail(r);
  }
  return res;
}

// The solution x of x B = v modulo p (B with independent rows), if any.
inline std::optional<ModVector> SolveRowsModP(ModMatrix const &B,
                                              ModVector const &v,
                                              T_residue p) {
  int r = B.rows();
  int c = B.cols();
  ModMatrix M(r + 1, c);
  M.topRows(r) = B;
  for (int j = 0; j < c; j++) {
    M(r, j) = v[j];
  }
  ModMatrix ker = LeftNullspaceModP(M, p);
  for (int u = 0; u < ker.rows(); u++) {
    if (ker(u, r) != 0) {
      T_residue inv = ResInv(ResNeg(ker(u, r), p), p);
      ModVector sol(r);
      for (int i = 0; i < r; i++) {
        sol[i] = ResMul(ker(u, i), inv, p);
      }
      return sol;
    }
  }
  return {};
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXMODPRIME_H_
// clang-format on
