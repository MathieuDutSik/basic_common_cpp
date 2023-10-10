// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_FUNCTIONS_H_
#define SRC_MATRIX_MAT_FUNCTIONS_H_

// clang-format off
#include "MAT_Matrix.h"
#include <set>
// clang-format on

//
// This routine is adequate for simplifying systems used
// for linear programming.
//
template <typename T> MyMatrix<T> SortUnicizeMatrix(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  auto comp = [](MyVector<T> const &V1, MyVector<T> const &V2) -> bool {
    return IsLower(V1, V2);
  };
  std::set<MyVector<T>, std::function<bool(MyVector<T>, MyVector<T>)>>
      eListVect(comp);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<T> V = GetMatrixRow(M, iRow);
    MyVector<T> Vcan = CanonicalizeVector(V);
    eListVect.insert(Vcan);
  }
  int nbVect = eListVect.size();
  MyMatrix<T> Mret(nbVect, nbCol);
  int idx = 0;
  for (auto &eV : eListVect) {
    for (int iCol = 0; iCol < nbCol; iCol++)
      Mret(idx, iCol) = eV(iCol);
    idx++;
  }
  return Mret;
}

// clang-format off
#endif  // SRC_MATRIX_MAT_FUNCTIONS_H_
// clang-format on
