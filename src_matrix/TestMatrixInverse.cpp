// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "QuadField.h"
#include "MAT_MatrixInt.h"
// clang-format on

template<typename T>
void process(int n) {
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
}


int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "TestMatrixInverse [arith] [n]\n";
      return -1;
    }
    std::string arith = argv[1];
    int n;
    sscanf(argv[2], "%d", &n);
    auto f=[&]() -> void {
      if (arith == "mpz_class")
        return process<mpz_class>(n);
      if (arith == "mpq_class")
        return process<mpq_class>(n);
      if (arith == "safe_integer")
        return process<SafeInt64>(n);
      if (arith == "safe_rational")
        return process<Rational<SafeInt64>>(n);
      std::cerr << "Failed to find a matching entry\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of TestMatrixInverse\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of TestMatrixInverse\n";
    exit(e.eVal);
  }
}
