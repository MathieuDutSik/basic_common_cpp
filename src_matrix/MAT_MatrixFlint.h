// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXFLINT_H_
#define SRC_MATRIX_MAT_MATRIXFLINT_H_

// The native flint matrix backend for MyMatrix<fmpz_class> and
// MyMatrix<fmpq_class>. fmpz_mat / fmpq_mat carry the fast algorithms
// (machine-word and multimodular products, modular determinants and
// Hermite / Smith normal forms with no coefficient explosion), so the
// generic entry points marshal into them, compute, and marshal back.
// The marshalling is O(n^2) single-word copies against the O(n^3)-and-up
// operation itself, measured below 3% of the cheapest routed operation.
//
// This header is only included when building with ENABLE_FLINT_SUPPORT.

// clang-format off
#include "NumberTheoryFlint.h"
#include "MatrixTypes.h"
#include <flint/fmpz_mat.h>
#include <flint/fmpq_mat.h>
#include <utility>
// clang-format on

inline void FlintSetMatrix(fmpz_mat_t out, MyMatrix<fmpz_class> const &M) {
  for (int i = 0; i < M.rows(); i++)
    for (int j = 0; j < M.cols(); j++)
      fmpz_set(fmpz_mat_entry(out, i, j), M(i, j).get_fmpz_t());
}

// The flint matrix is consumed: the entries are swapped out, so the heap
// limbs of the large values are handed over instead of copied.
inline MyMatrix<fmpz_class> FlintGetMatrix(fmpz_mat_t M) {
  MyMatrix<fmpz_class> ret(fmpz_mat_nrows(M), fmpz_mat_ncols(M));
  for (int i = 0; i < ret.rows(); i++)
    for (int j = 0; j < ret.cols(); j++)
      fmpz_swap(ret(i, j).get_fmpz_t(), fmpz_mat_entry(M, i, j));
  return ret;
}

inline void FlintSetMatrix(fmpq_mat_t out, MyMatrix<fmpq_class> const &M) {
  for (int i = 0; i < M.rows(); i++)
    for (int j = 0; j < M.cols(); j++)
      fmpq_set(fmpq_mat_entry(out, i, j), M(i, j).get_fmpq_t());
}

inline MyMatrix<fmpq_class> FlintGetMatrix(fmpq_mat_t M) {
  MyMatrix<fmpq_class> ret(fmpq_mat_nrows(M), fmpq_mat_ncols(M));
  for (int i = 0; i < ret.rows(); i++)
    for (int j = 0; j < ret.cols(); j++)
      fmpq_swap(ret(i, j).get_fmpq_t(), fmpq_mat_entry(M, i, j));
  return ret;
}

// Matrix product

inline MyMatrix<fmpz_class> FlintProductMatrix(MyMatrix<fmpz_class> const &A,
                                               MyMatrix<fmpz_class> const &B) {
  fmpz_mat_t a, b, c;
  fmpz_mat_init(a, A.rows(), A.cols());
  fmpz_mat_init(b, B.rows(), B.cols());
  fmpz_mat_init(c, A.rows(), B.cols());
  FlintSetMatrix(a, A);
  FlintSetMatrix(b, B);
  fmpz_mat_mul(c, a, b);
  MyMatrix<fmpz_class> C = FlintGetMatrix(c);
  fmpz_mat_clear(a);
  fmpz_mat_clear(b);
  fmpz_mat_clear(c);
  return C;
}

inline MyMatrix<fmpq_class> FlintProductMatrix(MyMatrix<fmpq_class> const &A,
                                               MyMatrix<fmpq_class> const &B) {
  fmpq_mat_t a, b, c;
  fmpq_mat_init(a, A.rows(), A.cols());
  fmpq_mat_init(b, B.rows(), B.cols());
  fmpq_mat_init(c, A.rows(), B.cols());
  FlintSetMatrix(a, A);
  FlintSetMatrix(b, B);
  fmpq_mat_mul(c, a, b);
  MyMatrix<fmpq_class> C = FlintGetMatrix(c);
  fmpq_mat_clear(a);
  fmpq_mat_clear(b);
  fmpq_mat_clear(c);
  return C;
}

// Determinant

inline fmpz_class FlintDeterminant(MyMatrix<fmpz_class> const &A) {
  fmpz_mat_t a;
  fmpz_mat_init(a, A.rows(), A.cols());
  FlintSetMatrix(a, A);
  fmpz_class det;
  fmpz_mat_det(det.get_fmpz_t(), a);
  fmpz_mat_clear(a);
  return det;
}

inline fmpq_class FlintDeterminant(MyMatrix<fmpq_class> const &A) {
  fmpq_mat_t a;
  fmpq_mat_init(a, A.rows(), A.cols());
  FlintSetMatrix(a, A);
  fmpq_class det;
  fmpq_mat_det(det.get_fmpq_t(), a);
  fmpq_mat_clear(a);
  return det;
}

// Row Hermite normal form, same convention as
// ComputeRowHermiteNormalForm_Kernel: H upper echelon with positive pivots,
// entries above a pivot in [0, pivot), zero rows at the bottom.

inline MyMatrix<fmpz_class>
FlintRowHermiteNormalForm(MyMatrix<fmpz_class> const &M) {
  fmpz_mat_t a, h;
  fmpz_mat_init(a, M.rows(), M.cols());
  fmpz_mat_init(h, M.rows(), M.cols());
  FlintSetMatrix(a, M);
  fmpz_mat_hnf(h, a);
  MyMatrix<fmpz_class> H = FlintGetMatrix(h);
  fmpz_mat_clear(a);
  fmpz_mat_clear(h);
  return H;
}

// Returns (U, H) with U unimodular and U M = H, matching the return of
// ComputeRowHermiteNormalForm. For a rank deficient M the transformation
// is not unique and flint may pick a different valid U than the generic
// kernel does; H is the same.
inline std::pair<MyMatrix<fmpz_class>, MyMatrix<fmpz_class>>
FlintRowHermiteNormalFormTransform(MyMatrix<fmpz_class> const &M) {
  fmpz_mat_t a, h, u;
  fmpz_mat_init(a, M.rows(), M.cols());
  fmpz_mat_init(h, M.rows(), M.cols());
  fmpz_mat_init(u, M.rows(), M.rows());
  FlintSetMatrix(a, M);
  fmpz_mat_hnf_transform(h, u, a);
  MyMatrix<fmpz_class> H = FlintGetMatrix(h);
  MyMatrix<fmpz_class> U = FlintGetMatrix(u);
  fmpz_mat_clear(a);
  fmpz_mat_clear(h);
  fmpz_mat_clear(u);
  return {std::move(U), std::move(H)};
}

// Column Hermite normal form through the transpose: the column convention
// of ComputeColHermiteNormalForm_Kernel is the exact mirror of the row one
// (M U = H, H lower echelon, positive pivots, entries left of a pivot in
// [0, pivot)), so H_col(M) = H_row(M^T)^T and U_col = U_row^T.

inline MyMatrix<fmpz_class>
FlintColHermiteNormalForm(MyMatrix<fmpz_class> const &M) {
  fmpz_mat_t a, at, h;
  fmpz_mat_init(a, M.rows(), M.cols());
  fmpz_mat_init(at, M.cols(), M.rows());
  fmpz_mat_init(h, M.cols(), M.rows());
  FlintSetMatrix(a, M);
  fmpz_mat_transpose(at, a);
  fmpz_mat_hnf(h, at);
  fmpz_mat_clear(a);
  fmpz_mat_clear(at);
  MyMatrix<fmpz_class> H(M.rows(), M.cols());
  for (int i = 0; i < H.rows(); i++)
    for (int j = 0; j < H.cols(); j++)
      fmpz_swap(H(i, j).get_fmpz_t(), fmpz_mat_entry(h, j, i));
  fmpz_mat_clear(h);
  return H;
}

// Returns (U, H) with U unimodular and M U = H, matching the return of
// ComputeColHermiteNormalForm.
inline std::pair<MyMatrix<fmpz_class>, MyMatrix<fmpz_class>>
FlintColHermiteNormalFormTransform(MyMatrix<fmpz_class> const &M) {
  fmpz_mat_t at, h, u;
  fmpz_mat_init(at, M.cols(), M.rows());
  fmpz_mat_init(h, M.cols(), M.rows());
  fmpz_mat_init(u, M.cols(), M.cols());
  {
    fmpz_mat_t a;
    fmpz_mat_init(a, M.rows(), M.cols());
    FlintSetMatrix(a, M);
    fmpz_mat_transpose(at, a);
    fmpz_mat_clear(a);
  }
  fmpz_mat_hnf_transform(h, u, at);
  MyMatrix<fmpz_class> H(M.rows(), M.cols());
  for (int i = 0; i < H.rows(); i++)
    for (int j = 0; j < H.cols(); j++)
      fmpz_swap(H(i, j).get_fmpz_t(), fmpz_mat_entry(h, j, i));
  MyMatrix<fmpz_class> U(M.cols(), M.cols());
  for (int i = 0; i < U.rows(); i++)
    for (int j = 0; j < U.cols(); j++)
      fmpz_swap(U(i, j).get_fmpz_t(), fmpz_mat_entry(u, j, i));
  fmpz_mat_clear(at);
  fmpz_mat_clear(h);
  fmpz_mat_clear(u);
  return {std::move(U), std::move(H)};
}

// Inverse. The rational inverse is unique, so the result is identical to
// the one of the generic kernel. For the integer type the inverse has to
// be representable over the integers, with the same failure behavior as
// InverseFractionFreeLU.

inline MyMatrix<fmpq_class> FlintInverse(MyMatrix<fmpq_class> const &A) {
  fmpq_mat_t a, b;
  fmpq_mat_init(a, A.rows(), A.cols());
  fmpq_mat_init(b, A.rows(), A.cols());
  FlintSetMatrix(a, A);
  int success = fmpq_mat_inv(b, a);
  fmpq_mat_clear(a);
  if (!success) {
    fmpq_mat_clear(b);
    std::cerr << "FlintInverse: the matrix is not invertible\n";
    throw TerminalException{1};
  }
  MyMatrix<fmpq_class> B = FlintGetMatrix(b);
  fmpq_mat_clear(b);
  return B;
}

inline MyMatrix<fmpz_class> FlintInverse(MyMatrix<fmpz_class> const &A) {
  fmpz_mat_t a, b;
  fmpz_mat_init(a, A.rows(), A.cols());
  fmpz_mat_init(b, A.rows(), A.cols());
  FlintSetMatrix(a, A);
  fmpz_t den;
  fmpz_init(den);
  int success = fmpz_mat_inv(b, den, a);
  fmpz_mat_clear(a);
  if (!success) {
    fmpz_mat_clear(b);
    fmpz_clear(den);
    std::cerr << "FlintInverse: the matrix is not invertible\n";
    throw TerminalException{1};
  }
  MyMatrix<fmpz_class> B(A.rows(), A.cols());
  for (int i = 0; i < B.rows(); i++)
    for (int j = 0; j < B.cols(); j++) {
      if (!fmpz_divisible(fmpz_mat_entry(b, i, j), den)) {
        fmpz_mat_clear(b);
        fmpz_clear(den);
        std::cerr << "FlintInverse: A^{-1} is not representable over the "
                     "integers\n";
        throw TerminalException{1};
      }
      fmpz_divexact(B(i, j).get_fmpz_t(), fmpz_mat_entry(b, i, j), den);
    }
  fmpz_mat_clear(b);
  fmpz_clear(den);
  return B;
}

// The kernel basis in the convention of NullspaceTrMat_Kernel: the reduced
// row echelon form is unique, and for every free column j the emitted
// vector has -1 at position j and the rref entry at each pivot column, so
// the output is identical to the one of the generic kernel.
inline MyMatrix<fmpq_class>
FlintNullspaceTrMat(MyMatrix<fmpq_class> const &A) {
  int nbRow = A.rows();
  int nbCol = A.cols();
  fmpq_mat_t a, r;
  fmpq_mat_init(a, nbRow, nbCol);
  fmpq_mat_init(r, nbRow, nbCol);
  FlintSetMatrix(a, A);
  slong rank = fmpq_mat_rref(r, a);
  fmpq_mat_clear(a);
  std::vector<int> ListColSelect;
  std::vector<uint8_t> ListColSelect01(nbCol, 0);
  for (slong i = 0; i < rank; i++)
    for (int j = 0; j < nbCol; j++)
      if (!fmpq_is_zero(fmpq_mat_entry(r, i, j))) {
        ListColSelect.push_back(j);
        ListColSelect01[j] = 1;
        break;
      }
  MyMatrix<fmpq_class> NSP =
      MyMatrix<fmpq_class>::Zero(nbCol - rank, nbCol);
  int nbVect = 0;
  for (int iCol = 0; iCol < nbCol; iCol++)
    if (ListColSelect01[iCol] == 0) {
      NSP(nbVect, iCol) = -1;
      for (slong iRank = 0; iRank < rank; iRank++)
        fmpq_set(NSP(nbVect, ListColSelect[iRank]).get_fmpq_t(),
                 fmpq_mat_entry(r, iRank, iCol));
      nbVect++;
    }
  fmpq_mat_clear(r);
  return NSP;
}

// The diagonal of the Smith normal form.

inline MyVector<fmpz_class>
FlintSmithNormalFormInvariant(MyMatrix<fmpz_class> const &M) {
  fmpz_mat_t a, s;
  fmpz_mat_init(a, M.rows(), M.cols());
  fmpz_mat_init(s, M.rows(), M.cols());
  FlintSetMatrix(a, M);
  fmpz_mat_snf(s, a);
  int siz = std::min(M.rows(), M.cols());
  MyVector<fmpz_class> V(siz);
  for (int i = 0; i < siz; i++)
    fmpz_swap(V(i).get_fmpz_t(), fmpz_mat_entry(s, i, i));
  fmpz_mat_clear(a);
  fmpz_mat_clear(s);
  return V;
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXFLINT_H_
// clang-format on
