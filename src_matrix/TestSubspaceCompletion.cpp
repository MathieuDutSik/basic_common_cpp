// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
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
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened\n";
    exit(e.eVal);
  }
  runtime(time);
}
