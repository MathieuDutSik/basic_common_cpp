// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format on

int main(int argc, char *argv[]) {
  using T = int;
  int n = 5;
  int m = 7;
  MyMatrix<T> eM(n, m);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++)
      eM(i, j) = 11 * i + j;
  T *ptr = eM.data();
  for (int i = 0; i < n; i++)
    for (int j = 0; j < m; j++) {
      T val1 = eM(i, j);
      int pos = j * n + i;
      T val2 = ptr[pos];
      if (val1 != val2) {
        std::cerr << "Inconsistency at i=" << i << " j=" << j << "\n";
      }
    }
}
