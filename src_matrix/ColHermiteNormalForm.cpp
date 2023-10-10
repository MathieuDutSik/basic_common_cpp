// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "MAT_MatrixInt.h"
// clang-format on

template <typename T>
void process(std::string const &FileI, std::string const &FileO) {
  MyMatrix<T> TheMat = ReadMatrixFile<T>(FileI);
  std::pair<MyMatrix<T>, MyMatrix<T>> ePair =
      ComputeColHermiteNormalForm(TheMat);
  //
  std::ofstream os(FileO);
  os << "return [";
  WriteMatrixGAP(os, ePair.first);
  os << ",\n";
  WriteMatrixGAP(os, ePair.second);
  os << "];\n";
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 4) {
      std::cerr << "This program is used as\n";
      std::cerr << "ColHermiteNormalForm arith [inputMat] [output]\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string FileO = argv[3];
    auto f = [&]() -> void {
      if (arith == "mpz_class")
        return process<mpz_class>(FileI, FileO);
      if (arith == "mpq_class")
        return process<mpq_class>(FileI, FileO);
      if (arith == "int")
        return process<int>(FileI, FileO);
      if (arith == "safe_integer")
        return process<SafeInt64>(FileI, FileO);
      if (arith == "safe_rational")
        return process<Rational<SafeInt64>>(FileI, FileO);
      std::cerr << "Failed to match\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong\n";
    exit(e.eVal);
  }
}
