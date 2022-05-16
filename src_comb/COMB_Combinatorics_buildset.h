// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_COMB_COMB_COMBINATORICS_BUILDSET_H_
#define SRC_COMB_COMB_COMBINATORICS_BUILDSET_H_

#include "MAT_Matrix.h"
#include <vector>

MyMatrix<int> BuildSet(int const &n, int const &Nval) {
  int TotalNb = MyPow<int>(Nval, n);
  std::vector<int> eVect(n, 0);
  auto fUpdate = [&]() -> int {
    for (int iVar = 0; iVar < n; iVar++)
      if (eVect[iVar] < Nval - 1) {
        eVect[iVar]++;
        for (int jVar = 0; jVar < iVar; jVar++)
          eVect[jVar] = 0;
        return 0;
      }
    return -1;
  };
  MyMatrix<int> RetMatrix(TotalNb, n);
  int idx = 0;
  while (true) {
    for (int i = 0; i < n; i++)
      RetMatrix(idx, i) = eVect[i];
    idx++;
    int test = fUpdate();
    if (test == -1)
      break;
  }
  return RetMatrix;
}

int PositionBuildSet(int const &n, int const &Nval, MyVector<int> const &V) {
  int pos = 0;
  int eProd = 1;
  for (int i = 0; i < n; i++) {
    pos += eProd * V(i);
    eProd *= Nval;
  }
  return pos;
}

// clang-format off
#endif  // SRC_COMB_COMB_COMBINATORICS_BUILDSET_H_
// clang-format on
