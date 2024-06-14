// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXFLT_H_
#define SRC_MATRIX_MAT_MATRIXFLT_H_

// clang-format off
#include <stdio.h>
#include <stdlib.h>
#include <vector>
// clang-format on

template <typename T>
std::vector<std::vector<T>>
InverseSquareMatrix(std::vector<std::vector<T>> TheMat) {
  int n, i, j, k, jSel;
  T MaxVal, hVal, eVal1, eVal2, eVal;
  std::vector<std::vector<T>> WorkMat, InvMat;
  std::vector<T> eVect;
  WorkMat = TheMat;
  n = TheMat.size();
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      if (i == j)
        eVal = 1;
      else
        eVal = 0;
      eVect.push_back(eVal);
    }
    InvMat.push_back(eVect);
    eVect.clear();
  }
  for (i = 0; i < n; i++) {
    MaxVal = abs(WorkMat[i][i]);
    jSel = i;
    for (j = i + 1; j < n; j++) {
      hVal = abs(WorkMat[i][j]);
      if (hVal > MaxVal) {
        jSel = j;
        MaxVal = hVal;
      }
    }
    if (i != jSel) {
      for (k = 0; k < n; k++) {
        eVal1 = WorkMat[k][i];
        eVal2 = WorkMat[k][jSel];
        WorkMat[k][i] = eVal1;
        WorkMat[k][jSel] = eVal2;
        eVal1 = InvMat[k][i];
        eVal2 = InvMat[k][jSel];
        InvMat[k][i] = eVal1;
        InvMat[k][jSel] = eVal2;
      }
    }
    T alpha = 1 / WorkMat[i][i];
    for (k = 0; k < n; k++) {
      WorkMat[k][i] = alpha * WorkMat[k][i];
      InvMat[k][i] = alpha * InvMat[k][i];
    }
    for (j = 0; j < n; j++)
      if (i != j) {
        T alpha2 = WorkMat[i][j];
        for (k = 0; k < n; k++) {
          WorkMat[k][j] = WorkMat[k][j] - WorkMat[k][i] * alpha2;
          InvMat[k][j] = InvMat[k][j] - InvMat[k][i] * alpha2;
        }
      }
  }
  return InvMat;
}

template <typename T>
void PrintEigenvalues(std::ostream &os, MyMatrix<T> const &eMat) {
  int n = eMat.rows();
  MyVector<T> ListEigVal(n);
  MyMatrix<T> ListEigVect(n, n);
  jacobi_double(eMat, ListEigVal, ListEigVect);
  for (int i = 0; i < n; i++) {
    T eEig = ListEigVal(i);
    os << "i=" << i << "/" << n << " eig=" << eEig << "\n";
  }
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXFLT_H_
// clang-format on
