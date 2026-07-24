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

template <typename T> MyMatrix<T> CongrMap(MyMatrix<T> const &eMat) {
  MyMatrix<T> TheInv = Inverse(eMat);
  return TransposedMat(TheInv);
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXINVERSE_H_
// clang-format on
