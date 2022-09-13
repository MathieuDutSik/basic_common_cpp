// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format on
int main(int argc, char *argv[]) {
  //  using T=mpq_class;
  using T = mpz_class;
  //  using T=int;
  //  using T=long;
  try {
    if (argc != 3) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "ColHermiteNormalForm [inputMat] [output]\n");
      return -1;
    }
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> TheMat = ReadMatrix<T>(INmat);
    std::pair<MyMatrix<T>, MyMatrix<T>> ePair =
        ComputeColHermiteNormalForm(TheMat);
    //
    std::ofstream os(argv[2]);
    os << "return [";
    WriteMatrixGAP(os, ePair.first);
    os << ",\n";
    WriteMatrixGAP(os, ePair.second);
    os << "];\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
