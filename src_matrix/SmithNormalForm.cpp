// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format on

template<typename T>
void process(std::string const& FileI) {
  MyMatrix<T> M = ReadMatrixFile<T>(FileI);
  // computing the Smith normal form
  std::pair<MyMatrix<T>, MyMatrix<T>> pair = SmithNormalForm(M);
  std::cerr << "ROW=\n";
  WriteMatrix(std::cerr, pair.first);
  std::cerr << "COL=\n";
  WriteMatrix(std::cerr, pair.second);
  //
  MyMatrix<T> RedMat = pair.first * M * pair.second;
  std::cerr << "RedMat=\n";
  WriteMatrix(std::cerr, RedMat);
}


int main(int argc, char *argv[]) {
  try {
    if (argc != 2) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "SmithNormalForm [arith] [inputMat]\n");
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    auto f=[&]() -> void {
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
