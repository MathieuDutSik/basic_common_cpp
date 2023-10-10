// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "MAT_MatrixInt.h"
// clang-format on

template <typename T> void process(std::string const &FileI, std::ostream &os) {
  MyVector<T> V1 = ReadVectorFile<T>(FileI);
  MyVector<T> V2 = ScalarCanonicalizationVector(V1);
  WriteVector(os, V2);
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "ScalarCanonicalization [arith] [inputVec]\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    auto f = [&](std::ostream &os) -> void {
      if (arith == "mpq_class") {
        using T = mpq_class;
        return process<T>(FileI, os);
      }
      if (arith == "mpz_class") {
        using T = mpz_class;
        return process<T>(FileI, os);
      }
      if (arith == "cpp_int") {
        using T = boost::multiprecision::cpp_int;
        return process<T>(FileI, os);
      }
      if (arith == "safe_integer") {
        using T = SafeInt64;
        return process<T>(FileI, os);
      }
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(FileI, os);
      }
      std::cerr << "Failed to find a matching type\n";
      throw TerminalException{1};
    };
    f(std::cerr);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
}
