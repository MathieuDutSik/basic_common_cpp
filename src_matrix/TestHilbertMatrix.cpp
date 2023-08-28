// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "rational.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "MAT_Matrix.h"
// clang-format on

template<typename T>
void process(int n) {
  //
  MyMatrix<T> eMat(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      T val(i + j + 1);
      eMat(i, j) = 1 / val;
    }
  }
  std::cerr << "eMat=\n";
  WriteMatrix(std::cerr, eMat);
  //
  MyMatrix<T> eInv = Inverse(eMat);
  std::cerr << "eInv=\n";
  WriteMatrix(std::cerr, eInv);
}



int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "TestHilbertMatrix is used as\n";
      std::cerr << "TestHilbertMatrix [arith] [n]\n";
      return -1;
    }
    std::string arith = argv[1];
    int n;
    sscanf(argv[2], "%d", &n);
    auto f=[&]() -> void {
      if (arith == "mpq_class")
        return process<mpq_class>(n);
      if (arith == "safe_rational")
        return process<Rational<SafeInt64>>(n);
      std::cerr << "Failed to find a matching type\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
}
