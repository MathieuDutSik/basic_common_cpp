// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXPROJECTION_H_
#define SRC_MATRIX_MAT_MATRIXPROJECTION_H_

// Determinant-based projection utilities: the orthogonal hyperplane of a set of
// vectors and the projector onto a subspace for a non-degenerate quadratic form.

// clang-format off
#include "MAT_MatrixDeterminant.h"
// clang-format on

template <typename T> MyVector<T> OrthogonalHyperplane(MyMatrix<T> const &M) {
  int nbVert = M.rows();
  MyVector<T> eVect(nbVert);
  MyMatrix<T> Mred(nbVert - 1, nbVert - 1);
  int eCoeff = 1;
  for (int iVert = 0; iVert < nbVert; iVert++) {
    int iRow = 0;
    for (int iLine = 0; iLine < nbVert; iLine++)
      if (iLine != iVert) {
        for (int iCol = 0; iCol < nbVert - 1; iCol++)
          Mred(iRow, iCol) = M(iLine, iCol);
        iRow++;
      }
    eVect(iVert) = eCoeff * DeterminantMat(Mred);
    eCoeff *= -1;
  }
  return eVect;
}

/*
  G has to be non-degenerate so that we can define the projector.
  --- "G"         is a (n x n) matrix.
  --- "Basis"     is a (p x n) basis of a matrix.
  --- "Basis * G" is a (p x n) matrix representing the (x,u_i) scalar products
  --- "Basis * G * Basis^T" is a matrix of the scalar products
  --- The final matrix is likely to be
      Basis^T (Basis * G * Basis^T) ^ (-1)  Basis G
  If p = n the formula simplifies to Identitity so the formula ought to be
  correct.
 */
template <typename T>
MyMatrix<T> GetProjectionMatrix(MyMatrix<T> const &G,
                                MyMatrix<T> const &Basis) {
  MyMatrix<T> Gred = Basis * G * Basis.transpose();
  if (DeterminantMat(Gred) == 0) {
    std::cerr << "G=\n";
    WriteMatrix(std::cerr, G);
    std::cerr << "Basis=\n";
    WriteMatrix(std::cerr, Basis);
    std::cerr << "Gred=\n";
    WriteMatrix(std::cerr, Gred);
    std::cerr << "The matrix Gred should be invertible\n";
    throw TerminalException{1};
  }
  MyMatrix<T> RetMat = Basis.transpose() * Inverse(Gred) * Basis * G;
  return RetMat;
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXPROJECTION_H_
// clang-format on
