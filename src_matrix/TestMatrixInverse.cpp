// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format on

int main(int argc, char *argv[]) {
  using T = mpq_class;
  //  using T=mpz_class;
  //  using T=int;
  //  using T=long;
  try {
    if (argc != 2) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      return -1;
    }
    // reading the matrix
    int n;
    sscanf(argv[1], "%d", &n);
    int nb = 100;
    for (int i = 0; i < nb; i++) {
      std::cerr << "i=" << i << "/" << nb << "\n";
      MyMatrix<T> eMat = RandomUnimodularMatrix<T>(n);
      MyMatrix<T> eInv = Inverse(eMat);
      MyMatrix<T> eProd = eMat * eInv;
      MyMatrix<T> eId = IdentityMat<T>(n);
      if (!TestEqualityMatrix(eProd, eId)) {
        std::cerr << "Error in matrix inverse\n";
        std::cerr << "eMat=\n";
        WriteMatrix(std::cerr, eMat);
        std::cerr << "eInv=\n";
        WriteMatrix(std::cerr, eInv);
        std::cerr << "eProd=\n";
        WriteMatrix(std::cerr, eProd);
        throw TerminalException{1};
      }
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
