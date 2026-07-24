// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXINVERSE_H_
#define SRC_MATRIX_MAT_MATRIXINVERSE_H_

// Matrix inverse over fields and rings.
//
//  --- TMat_Inverse_destroy: in-place Gauss-Jordan core (destroys its input).
//  --- InverseKernel / Inverse_destroy: thin wrappers over the core.
//  --- Inverse: dispatch entry point (field directly, other rings via the
//      overlying field).
//  --- CongrMap: the congruence map M -> (M^-1)^T.

// clang-format off
#include "MAT_MatrixFund.h"
// clang-format on

template <typename T>
void TMat_Inverse_destroy(MyMatrix<T> &Input, MyMatrix<T> &Output) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in TMat_Inverse_destroy");
  int iCol, iRow;
  int iRowB;
  int nbRow = Input.rows();
  int nbCol = Input.cols();
  T prov1;
  // Reuse-scratch for the row-combination product, hoisted so its buffer is
  // reused across the whole elimination. For fused-preferring types it collapses
  // to an empty object (no unused T is constructed); see is_fma_prefered.
  [[maybe_unused]] std::conditional_t<is_fma_prefered<T>::value, empty_scratch,
                                      T>
      prov2;
#ifdef DEBUG_MAT_MATRIX_DISABLE
  std::cerr << "TMat_Inverse_destroy, step 1\n";
#endif
#ifdef SANITY_CHECK_MAT_MATRIX
  if (nbRow != nbCol) {
    std::cerr << "Error on nbRow, nbCol in TMat_Inverse_destroy";
    throw TerminalException{1};
  }
#endif
  for (iCol = 0; iCol < nbRow; iCol++)
    for (iRow = 0; iRow < nbRow; iRow++) {
      if (iRow == iCol)
        prov1 = 1;
      else
        prov1 = 0;
      Output(iRow, iCol) = prov1;
    }
#ifdef DEBUG_MAT_MATRIX_DISABLE
  std::cerr << "TMat_Inverse_destroy, step 2\n";
#endif
  int iColFound;
  for (iRow = 0; iRow < nbRow; iRow++) {
#ifdef DEBUG_MAT_MATRIX_DISABLE
    std::cerr << "iRow=" << iRow << "\n";
    std::cerr << "Input=\n";
    WriteMatrix(std::cerr, Input);
#endif
    iColFound = SelectBestPivot<T>(
        nbCol, [&](int iCol) -> T const & { return Input(iRow, iCol); },
        [&](int iCol) -> bool { return iCol >= iRow; });
#ifdef SANITY_CHECK_MAT_MATRIX
    if (iColFound == -1) {
      std::cerr << "Error during the computation of the matrix inverse\n";
      throw TerminalException{1};
    }
#endif
    prov1 = 1 / Input(iRow, iColFound);
    for (iRowB = 0; iRowB < nbRow; iRowB++)
      Output(iRowB, iColFound) *= prov1;
    for (iRowB = iRow; iRowB < nbRow; iRowB++)
      Input(iRowB, iColFound) *= prov1;
    for (iCol = 0; iCol < nbCol; iCol++)
      if (iCol != iColFound) {
        prov1 = Input(iRow, iCol);
        if (prov1 != 0) {
          // Y -= prov1 * pivot_col. For types where the product allocates a
          // temporary each evaluation (mpz_class, mpq_class, ...) reuse prov2;
          // otherwise the fused form is at least as good (see is_fma_prefered).
          if constexpr (is_fma_prefered<T>::value) {
            for (iRowB = 0; iRowB < nbRow; iRowB++)
              Output(iRowB, iCol) -= prov1 * Output(iRowB, iColFound);
            for (iRowB = iRow; iRowB < nbRow; iRowB++)
              Input(iRowB, iCol) -= prov1 * Input(iRowB, iColFound);
          } else {
            for (iRowB = 0; iRowB < nbRow; iRowB++) {
              prov2 = prov1 * Output(iRowB, iColFound);
              Output(iRowB, iCol) -= prov2;
            }
            for (iRowB = iRow; iRowB < nbRow; iRowB++) {
              prov2 = prov1 * Input(iRowB, iColFound);
              Input(iRowB, iCol) -= prov2;
            }
          }
        }
      }
    if (iColFound != iRow) {
      for (iRowB = 0; iRowB < nbRow; iRowB++)
        std::swap(Output(iRowB, iColFound), Output(iRowB, iRow));
      for (iRowB = iRow; iRowB < nbRow; iRowB++)
        std::swap(Input(iRowB, iColFound), Input(iRowB, iRow));
    }
  }
#ifdef DEBUG_MAT_MATRIX_DISABLE
  std::cerr << "TMat_Inverse_destroy, step 3\n";
#endif
}

template <typename T> MyMatrix<T> InverseKernel(MyMatrix<T> const &Input) {
  int nbRow = Input.rows();
  MyMatrix<T> provMat = Input;
  MyMatrix<T> Output(nbRow, nbRow);
  TMat_Inverse_destroy(provMat, Output);
  return Output;
}

template <typename T> MyMatrix<T> Inverse_destroy(MyMatrix<T> &Input) {
  int nbRow = Input.rows();
  MyMatrix<T> Output(nbRow, nbRow);
  TMat_Inverse_destroy(Input, Output);
  return Output;
}

template <typename T>
requires is_ring_field<T>::value
inline MyMatrix<T> Inverse(MyMatrix<T> const &Input) {
  return InverseKernel(Input);
}

template <typename T>
requires (!is_ring_field<T>::value)
inline MyMatrix<T> Inverse(MyMatrix<T> const &Input) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> InputF = UniversalMatrixConversion<Tfield, T>(Input);
  MyMatrix<Tfield> OutputF = InverseKernel(InputF);
  return UniversalMatrixConversion<T, Tfield>(OutputF);
}

// Fraction-free Gauss-Jordan elimination, also known as fraction-free
// Gauss-Jordan or the Bareiss-Montante algorithm. Gauss-Jordan is run on the
// augmented matrix [A | I] using only ring operations and EXACT divisions: at
// step k every entry (pivot * M(i,j) - M(i,k) * M(k,j)) is divisible by the
// previous pivot, because Bareiss's theorem makes each intermediate entry a
// minor determinant of the input. This drives [A | I] to [det(A) * I | adj(A)],
// where adj(A) is the adjugate, so A^{-1} = adj(A) / det(A). Keeping the
// elimination fraction-free bounds intermediate operand growth (Hadamard), so
// for exact heavy arithmetic (rationals, number fields) it is faster than the
// classical Gauss-Jordan of InverseKernel, in the same way that
// DeterminantMatBareiss beats plain Gaussian elimination.
//
// Valid over any integral domain. The final division adj(A) / det(A) is exact
// over a field; over a ring of integers it is exact only when A is unimodular,
// and a non-exact division is reported rather than silently truncated.
template <typename T> MyMatrix<T> InverseBareiss(MyMatrix<T> const &Input) {
  int n = Input.rows();
  // Augmented matrix M = [A | I], of size n x 2n.
  MyMatrix<T> M(n, 2 * n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      M(i, j) = Input(i, j);
      M(i, n + j) = (i == j) ? T(1) : T(0);
    }
  T prev(1);
  for (int k = 0; k < n; k++) {
    if (M(k, k) == 0) {
      int r = -1;
      for (int i = k + 1; i < n; i++)
        if (M(i, k) != 0) {
          r = i;
          break;
        }
      if (r == -1) {
        std::cerr << "InverseBareiss: the matrix is singular\n";
        throw TerminalException{1};
      }
      M.row(k).swap(M.row(r));
    }
    T pivot = M(k, k);
    for (int i = 0; i < n; i++) {
      if (i == k)
        continue;
      for (int j = 0; j < 2 * n; j++) {
        if (j == k)
          continue;
        T val = pivot * M(i, j) - M(i, k) * M(k, j);
        T quot = val / prev; // exact division guaranteed by Bareiss's theorem
#ifdef DEBUG_MAT_MATRIX
        if (quot * prev != val) {
          std::cerr << "InverseBareiss: non-exact division, T is not an "
                       "integral domain\n";
          throw TerminalException{1};
        }
#endif
        M(i, j) = quot;
      }
      M(i, k) = T(0);
    }
    prev = pivot;
  }
  // Left block is now det(A) * I (the shared diagonal value) and the right
  // block is the adjugate adj(A) = det(A) * A^{-1}; the swap sign cancels in the
  // ratio, so no sign bookkeeping is needed.
  T det = prev;
  MyMatrix<T> Output(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      T num = M(i, n + j);
      T quot = num / det;
      if (quot * det != num) {
        std::cerr << "InverseBareiss: A^{-1} is not representable over T "
                     "(the determinant does not divide the adjugate)\n";
        throw TerminalException{1};
      }
      Output(i, j) = quot;
    }
  return Output;
}

template <typename T> MyMatrix<T> CongrMap(MyMatrix<T> const &eMat) {
  MyMatrix<T> TheInv = Inverse(eMat);
  return TransposedMat(TheInv);
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXINVERSE_H_
// clang-format on
