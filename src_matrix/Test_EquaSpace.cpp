// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format off

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "TestEqiaSpace [Space1] [Space2]\n";
      return -1;
    }
    // using T=mpz_class;
    using T = mpq_class;
    // reading the matrix
    std::string File1 = argv[1];
    std::string File2 = argv[2];
    MyMatrix<T> M1 = ReadMatrixFile<T>(File1);
    MyMatrix<T> M2 = ReadMatrixFile<T>(File2);
    bool test = TestEqualitySpannedSpaces(M1, M2);
    std::cerr << "test=" << test << "\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
