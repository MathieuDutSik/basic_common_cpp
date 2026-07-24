// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXDETERMINANT_H_
#define SRC_MATRIX_MAT_MATRIXDETERMINANT_H_

// Determinant computation over fields and rings.
//
// Several methods are collected here:
//  --- DeterminantMatKernel: Gaussian elimination over a field.
//  --- DeterminantMatBerkowitz: division-free Samuelson-Berkowitz, valid over
//      any commutative ring (including rings with zero divisors).
//  --- DeterminantMatUnitReduce: unit-pivot elimination plus a Berkowitz
//      residual, for rings opting into a division-free determinant (jets).
//  --- DeterminantMat: the dispatch entry point, selecting the right method
//      from the type traits of T.
//  --- DeterminantMatPermutation: the naive permutation-sum determinant, kept
//      as a slow control for consistency checks.

// clang-format off
#include "MAT_MatrixFund.h"
#include <algorithm>
#include <utility>
#include <vector>
// clang-format on

template <typename T> T DeterminantMatKernel(MyMatrix<T> const &TheMat) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in DeterminantMatKernel");
  T alpha;
  // Reuse-scratch for the elimination product (see is_fma_prefered); collapses
  // to an empty object for fused-preferring types.
  [[maybe_unused]]
  std::conditional_t<is_fma_prefered<T>::value, empty_scratch, T> prod;
  int n = TheMat.rows();
  MyMatrix<T> WorkMat = TheMat;
  std::vector<int> eVectPos(n, -1);
  T TheDet(1);
  for (int i = 0; i < n; i++) {
    int jSel = SelectBestPivot<T>(
        n, [&](int j) -> T const & { return WorkMat(i, j); },
        [&](int j) -> bool { return eVectPos[j] == -1; });
    if (jSel == -1) {
      return T(0);
    }
    eVectPos[jSel] = i;
    for (int j = 0; j < n; j++)
      if (j != jSel) {
        alpha = WorkMat(i, j) / WorkMat(i, jSel);
        if constexpr (is_fma_prefered<T>::value) {
          for (int k = 0; k < n; k++)
            WorkMat(k, j) -= alpha * WorkMat(k, jSel);
        } else {
          for (int k = 0; k < n; k++) {
            prod = alpha * WorkMat(k, jSel);
            WorkMat(k, j) -= prod;
          }
        }
      }
    TheDet = TheDet * WorkMat(i, jSel);
  }
  int nbchg = 0;
  for (int i = 0; i < n - 1; i++)
    for (int j = i + 1; j < n; j++)
      if (eVectPos[i] > eVectPos[j])
        nbchg++;
  int res = nbchg % 2;
  if (res == 0)
    return TheDet;
  return -TheDet;
}

// Bareiss fraction-free Gaussian elimination. Over an integral domain (e.g. the
// integers) it computes the determinant using only ring operations and EXACT
// divisions: at step k every entry M(i,j)*M(k,k) - M(i,k)*M(k,j) is divisible by
// the previous pivot, because Bareiss's theorem makes each intermediate entry a
// minor determinant of the input. So all quantities stay in the ring -- no
// fractions, and coefficient growth is bounded by the Hadamard bound rather than
// exploding -- while the cost stays O(n^3). This is the method of choice for
// exact integer matrices, where the field-based DeterminantMatKernel would drag
// in rational arithmetic. A zero pivot is handled by swapping in a non-zero entry
// from the same column below (flipping the sign); if none exists the matrix is
// singular and the determinant is zero.
template <typename T> T DeterminantMatBareiss(MyMatrix<T> const &Input) {
  int n = Input.rows();
  if (n == 0)
    return T(1);
  MyMatrix<T> M = Input;
  T prev(1);
  bool neg = false;
  for (int k = 0; k < n - 1; k++) {
    if (M(k, k) == 0) {
      int r = -1;
      for (int i = k + 1; i < n; i++)
        if (M(i, k) != 0) {
          r = i;
          break;
        }
      if (r == -1)
        return T(0);
      M.row(k).swap(M.row(r));
      neg = !neg;
    }
    for (int i = k + 1; i < n; i++)
      for (int j = k + 1; j < n; j++) {
        T val = M(i, j) * M(k, k) - M(i, k) * M(k, j);
        T quot = val / prev; // exact division guaranteed by Bareiss's theorem
#ifdef DEBUG_MAT_MATRIX
        if (quot * prev != val) {
          std::cerr << "DeterminantMatBareiss: non-exact division, T is not an "
                       "integral domain\n";
          throw TerminalException{1};
        }
#endif
        M(i, j) = quot;
      }
    prev = M(k, k);
  }
  T det = M(n - 1, n - 1);
  return neg ? -det : det;
}

// Samuelson-Berkowitz characteristic-polynomial determinant. It uses only +, -,
// * (no division), so it is valid over any commutative ring, including one with
// zero divisors such as the truncated jet ring. O(n^4). The characteristic
// polynomial coefficients p (p[0] = 1) are built as a product of lower-triangular
// Toeplitz matrices; det(A) = (-1)^n p[n].
template <typename T> T DeterminantMatBerkowitz(MyMatrix<T> const &A) {
  int n = A.rows();
  if (n == 0)
    return T(1);
  std::vector<T> p(1, T(1));
  for (int i = 1; i <= n; i++) {
    // Leading i x i submatrix partitioned as [[M, S], [R, a]] with
    // M = A[0..i-2][0..i-2], R = A(i-1, 0..i-2), S = A(0..i-2, i-1), a = A(i-1,i-1).
    // Toeplitz first column: c[0]=1, c[1]=-a, c[k+2] = -(R M^k S) for k=0..i-2.
    std::vector<T> c(i + 1);
    c[0] = T(1);
    c[1] = -A(i - 1, i - 1);
    if (i >= 2) {
      MyVector<T> y(i - 1); // y = M^k S, initialised to S
      for (int r = 0; r < i - 1; r++)
        y(r) = A(r, i - 1);
      for (int k = 0; k <= i - 2; k++) {
        T w(0); // w = R . y
        for (int r = 0; r < i - 1; r++)
          w += A(i - 1, r) * y(r);
        c[k + 2] = -w;
        if (k < i - 2) { // y <- M y
          MyVector<T> yn(i - 1);
          for (int r = 0; r < i - 1; r++) {
            T s(0);
            for (int col = 0; col < i - 1; col++)
              s += A(r, col) * y(col);
            yn(r) = s;
          }
          y = yn;
        }
      }
    }
    // new_p = (lower-triangular Toeplitz of c, (i+1) x i) * p
    std::vector<T> np(i + 1, T(0));
    for (int r = 0; r <= i; r++) {
      T s(0);
      for (int col = 0; col <= r && col < i; col++)
        s += c[r - col] * p[col];
      np[r] = s;
    }
    p = std::move(np);
  }
  T det = p[n];
  if (n % 2 == 1)
    det = -det;
  return det;
}

// Determinant over a ring with zero divisors (opting into
// determinant_division_free, e.g. jets). Eliminate with UNIT pivots -- entries
// whose constant term is non-zero, so the division is exact -- as far as
// possible, an ordinary O(n^3) Gaussian elimination. When the active block has no
// unit pivot left it is the corank-d residual (the matrix is singular at the
// degeneracy point), all of whose entries are zero divisors; its determinant is
// computed division-free by Berkowitz. Total O(n^3 + d^4): a non-singular matrix
// has no residual and this is a plain elimination; only the small corank-d block
// pays the division-free price.
template <typename T> T DeterminantMatUnitReduce(MyMatrix<T> const &Input) {
  using Tct = std::decay_t<decltype(constant_term(std::declval<T const &>()))>;
  int n = Input.rows();
  MyMatrix<T> M = Input;
  T det(1);
  bool neg = false;
  for (int step = 0; step < n; step++) {
    int pr = -1, pc = -1;
    decltype(f_cost_pivot(std::declval<T>())) best_cost{};
    for (int i = step; i < n; i++)
      for (int j = step; j < n; j++)
        if (constant_term(M(i, j)) != Tct(0)) { // unit (invertible) pivot
          auto cost = f_cost_pivot(M(i, j));
          if (pr == -1 || is_preferable_pivot(cost, best_cost)) {
            pr = i;
            pc = j;
            best_cost = cost;
          }
        }
    if (pr == -1) {
      int d = n - step;
      MyMatrix<T> R(d, d);
      for (int i = 0; i < d; i++)
        for (int j = 0; j < d; j++)
          R(i, j) = M(step + i, step + j);
      det = det * DeterminantMatBerkowitz(R);
      return neg ? -det : det;
    }
    if (pr != step) {
      M.row(step).swap(M.row(pr));
      neg = !neg;
    }
    if (pc != step) {
      M.col(step).swap(M.col(pc));
      neg = !neg;
    }
    det = det * M(step, step);
    for (int i = step + 1; i < n; i++) {
      T factor = M(i, step) / M(step, step); // unit pivot -> exact division
      for (int j = step + 1; j < n; j++)
        M(i, j) -= factor * M(step, j);
    }
  }
  return neg ? -det : det;
}

template <typename T>
requires (is_ring_field<T>::value && !determinant_division_free<T>::value)
inline T DeterminantMat(MyMatrix<T> const &Input) {
  return DeterminantMatKernel(Input);
}

// Determinant over a field that carries zero divisors and opts into a
// division-free algorithm (jets): elimination with unit pivots plus a Berkowitz
// determinant on the corank-d residual (see DeterminantMatUnitReduce).
template <typename T>
requires (is_ring_field<T>::value && determinant_division_free<T>::value)
inline T DeterminantMat(MyMatrix<T> const &Input) {
  return DeterminantMatUnitReduce(Input);
}

// Rings of integers (integral domains with exact division): compute the
// determinant with Bareiss fraction-free elimination, which stays inside the
// ring. This is dramatically faster than mapping to the fraction field and
// running rational Gaussian elimination -- benchmarks show more than 20x at
// n=200 -- because the rational path pays a GCD reduction at every elimination
// step whereas Bareiss keeps every entry integral and bounded by the minor
// (Hadamard) size.
template <typename T>
requires (!is_ring_field<T>::value && is_implementation_of_Z<T>::value)
inline T DeterminantMat(MyMatrix<T> const &Input) {
  return DeterminantMatBareiss(Input);
}

// Any other non-field ring: map to the overlying field and eliminate there.
template <typename T>
requires (!is_ring_field<T>::value && !is_implementation_of_Z<T>::value)
inline T DeterminantMat(MyMatrix<T> const &Input) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> InputF = UniversalMatrixConversion<Tfield, T>(Input);
  Tfield eDet_field = DeterminantMatKernel(InputF);
  return UniversalScalarConversion<T, Tfield>(eDet_field);
}

// A significantly slower algorithm for computing the determinant.
// It is good as a control for the above method and can also be used
// for consistency checks of arithmetics.
template <typename T> T DeterminantMatPermutation(MyMatrix<T> const &A) {
  int n = A.rows();
  if (n == 0)
    return T(1);
  std::vector<int> s(n);
  for (int i = 0; i < n; i++)
    s[i] = i;
  T TheDet(0);
  do {
    T eProd(1);
    for (int u = 0; u < n; u++)
      eProd *= A(u, s[u]);
    int eSign = 1;
    for (int i = 0; i < n; i++)
      for (int j = i + 1; j < n; j++)
        if (s[j] < s[i])
          eSign = -eSign;
    TheDet += eSign * eProd;
  } while (std::next_permutation(s.begin(), s.end()));
  return TheDet;
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXDETERMINANT_H_
// clang-format on
