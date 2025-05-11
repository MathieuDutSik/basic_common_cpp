// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "MAT_MatrixInt.h"
// clang-format on

template <typename T> void process(std::string const &FileI) {
  MyMatrix<T> M = ReadMatrixFile<T>(FileI);
  // computing the Smith normal form
  MyVector<T> VectInv = SmithNormalFormInvariant(M);
  //
  std::cerr << "VectInv=\n";
  WriteVector(std::cerr, VectInv);
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 2) {
      std::cerr << "This program is used as\n";
      std::cerr << "SmithNormalForm [arith] [inputMat]\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    auto f = [&]() -> void {
      if (arith == "mpz_class")
        return process<mpz_class>(FileI);
      if (arith == "mpq_class")
        return process<mpq_class>(FileI);
      if (arith == "safe_integer")
        return process<SafeInt64>(FileI);
      if (arith == "safe_rational")
        return process<Rational<SafeInt64>>(FileI);
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
