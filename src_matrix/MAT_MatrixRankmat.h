// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXRANKMAT_H_
#define SRC_MATRIX_MAT_MATRIXRANKMAT_H_

// The rank of a matrix: Gaussian elimination over a field, fraction-free
// Bareiss elimination over a ring.

// clang-format off
#include "MAT_MatrixFund.h"
// clang-format on

template <typename T> int RankMatKernel(MyMatrix<T> const &Input) {
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(Input);
  return eSelect.TheRank;
}

template <typename T>
requires is_ring_field<T>::value
inline int RankMat(MyMatrix<T> const &Input) {
  return RankMatKernel(Input);
}

// Fraction-free rank computation over a ring (Bareiss). The elimination
//     M(i,j) <- (pivot * M(i,j) - M(i,pivcol) * M(prow,j)) / prev
// divides by the previous pivot, which is exact over any integral
// domain (every entry stays a minor of the input), so no conversion to
// the overlying field and no gcd computations are needed, and the entry
// sizes are bounded by the minors instead of blowing up.
template <typename T>
requires (!is_ring_field<T>::value)
inline int RankMat(MyMatrix<T> const &Input) {
  int nbRow = Input.rows();
  int nbCol = Input.cols();
  MyMatrix<T> M = Input;
  T prev(1);
  int rank = 0;
  while (rank < nbRow && rank < nbCol) {
    // A non-zero pivot in the remaining block, swapped into position.
    int pRow = -1;
    int pCol = -1;
    for (int iRow = rank; iRow < nbRow && pRow < 0; iRow++) {
      for (int iCol = rank; iCol < nbCol && pRow < 0; iCol++) {
        if (M(iRow, iCol) != 0) {
          pRow = iRow;
          pCol = iCol;
        }
      }
    }
    if (pRow < 0) {
      break;
    }
    if (pRow != rank) {
      M.row(rank).swap(M.row(pRow));
    }
    if (pCol != rank) {
      M.col(rank).swap(M.col(pCol));
    }
    T pivot = M(rank, rank);
    for (int iRow = rank + 1; iRow < nbRow; iRow++) {
      for (int iCol = rank + 1; iCol < nbCol; iCol++) {
        T val = pivot * M(iRow, iCol) - M(iRow, rank) * M(rank, iCol);
        M(iRow, iCol) = val / prev; // exact division (Bareiss)
      }
      M(iRow, rank) = T(0);
    }
    prev = pivot;
    rank++;
  }
  return rank;
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXRANKMAT_H_
// clang-format on
