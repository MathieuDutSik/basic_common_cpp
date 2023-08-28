// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format off

template<typename T>
void process(std::string const& FileI) {
  // reading the matrix
  MyMatrix<T> M = ReadMatrixFile<T>(FileI);
  std::cerr << "M=\n";
  WriteMatrix(std::cerr, M);
  //
  FractionMatrix<T> eRec = RemoveFractionMatrixPlusCoeff(M);
  std::cerr << "Mult=" << eRec.TheMult << "\n";
  std::cerr << "TheMat=\n";
  WriteMatrix(std::cerr, eRec.TheMat);
}


int main(int argc, char *argv[]) {
  using T = mpq_class;
  try {
    if (argc != 2) {
      std::cerr << "This program is used as\n";
      std::cerr << "RemoveFractionMatrix [arith] [inputMat]\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    auto f=[&]() -> void {
      if (arith == "mpq_class")
        process<mpq_class>(FileI);
      if (arith == "safe_rational")
        process<Rational<SafeInt64>>(FileI);
      std::cerr << "Failed to find a matching type\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of RemoveFractionMatrix\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of RemoveFractionMatrix\n";
    exit(e.eVal);
  }
}
