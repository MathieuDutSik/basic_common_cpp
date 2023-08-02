// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format on


template<typename T>
void process(std::string const& FileI, std::string const& FileO) {
  MyMatrix<T> TheMat = ReadMatrixFile<T>(FileI);
  std::pair<MyMatrix<T>, MyMatrix<T>> ePair =
    ComputeRowHermiteNormalForm(TheMat);
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
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "RowHermiteNormalForm [arith] [inputMat] [output]\n");
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string FileO = argv[3];
    auto f=[&]() -> void {
      if (arith == "mpq_class")
        return process<mpq_class>(FileI, FileO);
      if (arith == "mpz_class")
        return process<mpz_class>(FileI, FileO);
      if (arith == "cpp_int")
        return process<cpp_int>(FileI, FileO);
      if (arith == "safe_integer")
        return process<SafeInt64>(FileI, FileO);
      if (arith == "safe_rational")
        return process<Rational<SafeInt64>>(FileI, FileO);
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
