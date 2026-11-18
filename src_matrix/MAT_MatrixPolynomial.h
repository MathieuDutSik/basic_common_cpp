// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXPOLYNOMIAL_H_
#define SRC_MATRIX_MAT_MATRIXPOLYNOMIAL_H_

// clang-format off
#include "MAT_Matrix.h"
#include "MAT_MatrixDeterminant.h"
#include "MAT_MatrixSolutionMat.h"
#include <optional>
#include <vector>
// clang-format on

/*
  Polynomials over a field T as the vectors of their coefficients, of
  increasing degree, and their action on square matrices.
 */

template <typename T>
std::vector<T> PolyProductCoeffs(std::vector<T> const &a,
                                 std::vector<T> const &b) {
  std::vector<T> c(a.size() + b.size() - 1, T(0));
  for (size_t i = 0; i < a.size(); i++) {
    for (size_t j = 0; j < b.size(); j++) {
      c[i + j] += a[i] * b[j];
    }
  }
  return c;
}

// The resultant of a and b, the determinant of their Sylvester matrix.
template <typename T>
T PolyResultant(std::vector<T> const &a, std::vector<T> const &b) {
  int da = static_cast<int>(a.size()) - 1;
  int db = static_cast<int>(b.size()) - 1;
  if (da == 0 && db == 0) {
    return T(1);
  }
  if (da == 0) {
    return MyPow<T>(a[0], db);
  }
  if (db == 0) {
    return MyPow<T>(b[0], da);
  }
  int m = da + db;
  MyMatrix<T> S = ZeroMatrix<T>(m, m);
  for (int i = 0; i < db; i++) {
    for (int k = 0; k <= da; k++) {
      S(i, i + k) = a[da - k];
    }
  }
  for (int i = 0; i < da; i++) {
    for (int k = 0; k <= db; k++) {
      S(db + i, i + k) = b[db - k];
    }
  }
  return DeterminantMat(S);
}

// M^e, e >= 0, by binary exponentiation.
template <typename T> MyMatrix<T> MatrixPower(MyMatrix<T> const &M, int64_t e) {
  MyMatrix<T> res = IdentityMat<T>(M.rows());
  MyMatrix<T> b = M;
  while (e > 0) {
    if (e & 1) {
      res = res * b;
    }
    e >>= 1;
    if (e > 0) {
      b = b * b;
    }
  }
  return res;
}

// p(M) for the polynomial p with integer coefficients, by Horner.
template <typename T>
MyMatrix<T> EvaluatePolynomialMatrix(std::vector<int64_t> const &coeff,
                                     MyMatrix<T> const &M) {
  int n = M.rows();
  MyMatrix<T> res = ZeroMatrix<T>(n, n);
  for (size_t k = coeff.size(); k-- > 0;) {
    res = res * M;
    T c = UniversalScalarConversion<T, int64_t>(coeff[k]);
    for (int i = 0; i < n; i++) {
      res(i, i) += c;
    }
  }
  return res;
}

// The monic minimal polynomial of the square matrix g: the first power of
// g in the span of the lower ones.
template <typename T>
std::vector<T> MinimalPolynomialCoeffs(MyMatrix<T> const &g) {
  int n = g.rows();
  auto flat = [&](MyMatrix<T> const &M) -> MyVector<T> {
    MyVector<T> v(n * n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        v(i * n + j) = M(i, j);
      }
    }
    return v;
  };
  std::vector<MyVector<T>> pw{flat(IdentityMat<T>(n))};
  MyMatrix<T> M = IdentityMat<T>(n);
  while (true) {
    M = M * g;
    MyVector<T> w = flat(M);
    std::optional<MyVector<T>> opt =
        SolutionMat(MatrixFromVectorFamily(pw), w);
    if (opt) {
      // g^k = sum_i c_i g^i
      std::vector<T> mu;
      for (int i = 0; i < opt->size(); i++) {
        mu.push_back(-(*opt)(i));
      }
      mu.push_back(T(1));
      return mu;
    }
    pw.push_back(w);
  }
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXPOLYNOMIAL_H_
// clang-format on
