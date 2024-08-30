// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "rational.h"
#include "MAT_Matrix.h"
// clang-format on

int main(int argc, char *argv[]) {
  //  using T=mpq_class;
  using T = Rational<int64_t>;
  try {
    if (argc != 3) {
      std::cerr << "TestProdMat [M1] [M2]\n";
      return -1;
    }
    // reading the matrix
    std::string file1 = argv[1];
    std::string file2 = argv[2];
    MyMatrix<T> M1 = ReadMatrixFile<T>(file1);
    MyMatrix<T> M2 = ReadMatrixFile<T>(file2);
    std::cerr << "M1=\n";
    WriteMatrix(std::cerr, M1);
    std::cerr << "M2=\n";
    WriteMatrix(std::cerr, M2);
    //
    MyMatrix<T> eProd = M1 * M2;
    std::cerr << "eProd=\n";
    WriteMatrix(std::cerr, eProd);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
