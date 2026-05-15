// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format on

template <typename T> void process(int n, int m) {
  int nb = 100;
  for (int i = 0; i < nb; i++) {
    MyMatrix<T> Munimodular = RandomUnimodularMatrix<T>(n);
    std::vector<MyVector<T>> l_vect;
    for (int i = 0; i < n - m; i++)
      l_vect.push_back(GetMatrixRow(Munimodular, i));
    MyMatrix<T> M = MatrixFromVectorFamily(l_vect);
    MyMatrix<T> Msub = SubspaceCompletionInt(M, n);
    MyMatrix<T> Mconcat = Concatenate(M, Msub);
    T eDet = DeterminantMat(Mconcat);
    if (T_abs(eDet) != 1) {
      std::cerr << "eDet = " << eDet << " while it should be 1 or -1\n";
      throw TerminalException{1};
    }
  }
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4) {
      std::cerr << "Test_SubspaceCompletion is used as\n";
      std::cerr << "Test_SubspaceCompletion [arith] [n] [m]\n";
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
    std::cerr << "Something wrong happened\n";
    exit(e.eVal);
  }
  runtime(time);
}
