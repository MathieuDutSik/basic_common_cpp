// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXMODINT64_H_
#define SRC_MATRIX_MAT_MATRIXMODINT64_H_

// clang-format off
#include "MAT_SparseMatrixMod.h"
#include <cstdint>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

/*
  Small dense matrices modulo m < 2^31 (so that the products of two
  residues fit in an int64_t), stored row major in a std::vector<int64_t>:
  no allocation per entry, for the inner loops (the coset keys of the
  congruence subgroups compute an inverse modulo N per key).

  For a prime m the elimination divides by the pivots, as
  DeterminantMatMod does. For a composite m a nonzero pivot need not be
  invertible, and the entries below the diagonal are cleared by the
  Euclidean algorithm on the rows instead, which uses no inverse: the
  diagonal entry of a column is then a gcd of its entries.
 */

// The elimination to an upper triangular matrix of the n x n matrix a,
// with the same row operations applied to b (which can be empty). Returns
// the sign (+-1) of the row permutations, or 0 if a column has no
// invertible pivot for a prime m (then det = 0).
inline int64_t EliminationMatModInt64(std::vector<int64_t> &a, std::vector<int64_t> &b, int const &n,
                                      int64_t const &m, bool const &is_prime) {
  bool has_b = !b.empty();
  int64_t sign = 1;
  auto f_swap = [&](int const &i1, int const &i2) -> void {
    for (int j = 0; j < n; j++) {
      std::swap(a[i1 * n + j], a[i2 * n + j]);
    }
    if (has_b) {
      for (int j = 0; j < n; j++) {
        std::swap(b[i1 * n + j], b[i2 * n + j]);
      }
    }
    sign = -sign;
  };
  // row i_dst -= q * row i_src, on a from the column col (the entries
  // before are zero in both rows) and on all of b.
  auto f_sub = [&](int const &i_dst, int64_t const &q, int const &i_src, int const &col) -> void {
    if (q == 0) {
      return;
    }
    int64_t mq = m - q;
    for (int j = col; j < n; j++) {
      a[i_dst * n + j] = (a[i_dst * n + j] + mq * a[i_src * n + j]) % m;
    }
    if (has_b) {
      for (int j = 0; j < n; j++) {
        b[i_dst * n + j] = (b[i_dst * n + j] + mq * b[i_src * n + j]) % m;
      }
    }
  };
  for (int col = 0; col < n; col++) {
    if (is_prime) {
      int piv = -1;
      for (int r = col; r < n; r++) {
        if (a[r * n + col] != 0) {
          piv = r;
          break;
        }
      }
      if (piv == -1) {
        return 0;
      }
      if (piv != col) {
        f_swap(piv, col);
      }
      int64_t inv = InverseModPrime(a[col * n + col], m);
      for (int r = col + 1; r < n; r++) {
        f_sub(r, (a[r * n + col] * inv) % m, col, col);
      }
    } else {
      for (int r = col + 1; r < n; r++) {
        while (a[r * n + col] != 0) {
          f_sub(col, a[col * n + col] / a[r * n + col], r, col);
          f_swap(col, r);
        }
      }
    }
  }
  return sign;
}

// The determinant modulo m of the n x n matrix a (row major, entries in
// [0, m)).
inline int64_t DeterminantMatModInt64(std::vector<int64_t> const &a, int const &n, int64_t const &m,
                                      bool const &is_prime) {
  if (m == 1) {
    return 0;
  }
  std::vector<int64_t> a_work = a;
  std::vector<int64_t> b_empty;
  int64_t sign = EliminationMatModInt64(a_work, b_empty, n, m, is_prime);
  if (sign == 0) {
    return 0;
  }
  int64_t det = ResidueModPrime(sign, m);
  for (int i = 0; i < n; i++) {
    det = (det * a_work[i * n + i]) % m;
  }
  return det;
}

// The inverse modulo m of the n x n matrix a (row major, entries in
// [0, m)), empty if it is not invertible.
inline std::optional<std::vector<int64_t>> InverseMatModInt64(std::vector<int64_t> const &a, int const &n,
                                                              int64_t const &m, bool const &is_prime) {
  std::vector<int64_t> inv(n * n, 0);
  if (m == 1) {
    return inv;
  }
  std::vector<int64_t> a_work = a;
  for (int i = 0; i < n; i++) {
    inv[i * n + i] = 1;
  }
  if (EliminationMatModInt64(a_work, inv, n, m, is_prime) == 0) {
    return {};
  }
  // The back substitution, from the last row.
  for (int col = n - 1; col >= 0; col--) {
    int64_t piv = a_work[col * n + col];
    if (std::gcd(piv, m) != 1) {
      return {};
    }
    int64_t piv_inv = InverseModPrime(piv, m);
    for (int j = col; j < n; j++) {
      a_work[col * n + j] = (a_work[col * n + j] * piv_inv) % m;
    }
    for (int j = 0; j < n; j++) {
      inv[col * n + j] = (inv[col * n + j] * piv_inv) % m;
    }
    for (int r = 0; r < col; r++) {
      int64_t q = a_work[r * n + col];
      if (q != 0) {
        int64_t mq = m - q;
        a_work[r * n + col] = 0;
        for (int j = 0; j < n; j++) {
          inv[r * n + j] = (inv[r * n + j] + mq * inv[col * n + j]) % m;
        }
      }
    }
  }
  return inv;
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXMODINT64_H_
// clang-format on
