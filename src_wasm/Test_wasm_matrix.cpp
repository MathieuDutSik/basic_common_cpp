// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm matrix test: exercises MyMatrix<cpp_rational> / MyMatrix<cpp_int>
// for inverse, determinant, nullspace, and integer-side HNF.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include "MAT_MatrixInt.h"
#include <iostream>

using cpp_int = boost::multiprecision::cpp_int;
using cpp_rational = boost::multiprecision::cpp_rational;

static int n_fail = 0;

#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << " : " << #cond    \
                << "\n";                                                       \
      n_fail++;                                                                \
    }                                                                          \
  } while (0)

// Inverse roundtrip: random unimodular matrix M, eInv = Inverse(M).
// Then M * eInv must equal the identity.
static void test_unimodular_inverse() {
  for (int trial = 0; trial < 30; trial++) {
    MyMatrix<cpp_int> M_int = RandomUnimodularMatrix<cpp_int>(4);
    MyMatrix<cpp_rational> M(4, 4);
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        M(i, j) = cpp_rational(M_int(i, j));
    MyMatrix<cpp_rational> eInv = Inverse(M);
    MyMatrix<cpp_rational> eProd = M * eInv;
    MyMatrix<cpp_rational> eId = IdentityMat<cpp_rational>(4);
    CHECK(TestEqualityMatrix(eProd, eId));
  }
}

// Determinant of an integer unimodular matrix is ±1.
static void test_determinant_unimodular() {
  for (int trial = 0; trial < 30; trial++) {
    MyMatrix<cpp_int> M = RandomUnimodularMatrix<cpp_int>(5);
    cpp_int det = DeterminantMat(M);
    CHECK(det == 1 || det == -1);
  }
}

// Hilbert-matrix style rational matrix: build H_ij = 1/(i+j+1) and confirm
// H * H^{-1} == I.
static void test_hilbert_inverse() {
  int n = 4;
  MyMatrix<cpp_rational> H(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      H(i, j) = cpp_rational(1, i + j + 1);
  MyMatrix<cpp_rational> Hinv = Inverse(H);
  MyMatrix<cpp_rational> Prod = H * Hinv;
  MyMatrix<cpp_rational> Id = IdentityMat<cpp_rational>(n);
  CHECK(TestEqualityMatrix(Prod, Id));
}

// Nullspace property: NullspaceMat(M) * M = 0  (rows of NullspaceMat span
// the left nullspace of M).
static void test_nullspace_property() {
  int n_row = 5, n_col = 3;
  MyMatrix<cpp_rational> M(n_row, n_col);
  // A rank-2 matrix whose left nullspace is 3-dimensional.
  M(0, 0) = 1; M(0, 1) = 0; M(0, 2) = 1;
  M(1, 0) = 0; M(1, 1) = 1; M(1, 2) = 2;
  M(2, 0) = 1; M(2, 1) = 1; M(2, 2) = 3;
  M(3, 0) = 2; M(3, 1) = -1; M(3, 2) = 0;
  M(4, 0) = -1; M(4, 1) = 3; M(4, 2) = 5;
  MyMatrix<cpp_rational> NS = NullspaceMat(M);
  MyMatrix<cpp_rational> Prod = NS * M;
  CHECK(Prod.rows() == NS.rows());
  CHECK(Prod.cols() == n_col);
  for (int i = 0; i < Prod.rows(); i++)
    for (int j = 0; j < Prod.cols(); j++)
      CHECK(Prod(i, j) == 0);
  CHECK(NS.rows() == n_row - 2);
}

// HNF roundtrip: ComputeRowHermiteNormalForm returns (P, H) with P*M = H.
static void test_hnf_roundtrip() {
  for (int trial = 0; trial < 10; trial++) {
    int n = 4, m = 4;
    MyMatrix<cpp_int> M(n, m);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < m; j++)
        M(i, j) = cpp_int(rand() % 11 - 5);
    auto [P, H] = ComputeRowHermiteNormalForm(M);
    MyMatrix<cpp_int> ProdH = P * M;
    CHECK(TestEqualityMatrix(ProdH, H));
    CHECK(T_abs(DeterminantMat(P)) == 1);
  }
}

int main() {
  test_unimodular_inverse();
  test_determinant_unimodular();
  test_hilbert_inverse();
  test_nullspace_property();
  test_hnf_roundtrip();
  if (n_fail == 0) {
    std::cerr << "Test_wasm_matrix: all checks passed.\n";
    return 0;
  }
  std::cerr << "Test_wasm_matrix: " << n_fail << " check(s) failed.\n";
  return 1;
}
