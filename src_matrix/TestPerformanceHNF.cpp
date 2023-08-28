// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format on

int main(int argc, char *argv[]) {
  //  using T=mpq_class;
  using T = mpz_class;
  //  using T=int;
  //  using T=long;
  try {
    if (argc != 3) {
      std::cerr << "TestPerformanceHNF is used as\n";
      std::cerr << "TestPerformanceHNF [n] [m]\n";
      return -1;
    }
    // reading the matrix
    int n, m;
    sscanf(argv[1], "%d", &n);
    sscanf(argv[2], "%d", &m);
    int nb = 100;
    int siz = 2;
    for (int i = 0; i < nb; i++) {
      std::cerr << "i=" << i << "/" << nb << "\n";
      MyMatrix<T> eMat(n, m);
      for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
          eMat(i, j) = random() % (2 * siz + 1) - siz;
      std::pair<MyMatrix<T>, MyMatrix<T>> ePair =
          ComputeRowHermiteNormalForm(eMat);
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
