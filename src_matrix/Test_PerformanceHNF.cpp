// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format on

template <typename T> void process(int n, int m) {
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
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4) {
      std::cerr << "TestPerformanceHNF is used as\n";
      std::cerr << "TestPerformanceHNF [arith] [n] [m]\n";
      return -1;
    }
    std::string arith = argv[1];
    int n, m;
    sscanf(argv[2], "%d", &n);
    sscanf(argv[3], "%d", &m);
    if (arith == "mpz_class") {
      process<mpz_class>(n, m);
    } else if (arith == "SafeInt64") {
      process<SafeInt64>(n, m);
    } else {
      std::cerr << "Unknown arith: " << arith << "\n";
      throw TerminalException{1};
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
