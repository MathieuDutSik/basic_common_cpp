// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "MAT_MatrixInt.h"
// clang-format on

template <typename T> void process(std::string const& FileI) {
  std::vector<MyMatrix<T>> ListMat = ReadListMatrixFile<T>(FileI);
  WriteListMatrix(std::cerr, ListMat);
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "Test_ReadListMat [arith] [FileI]\n";
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
      std::cerr << "Failed to find a matching entry for arith=" << arith << "\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of Test_ReadListMat\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Test_ReadListMat\n";
    exit(e.eVal);
  }
  runtime(time);
}
