// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format off

int main(int argc, char *argv[]) {
  HumanTime time;
  using T = mpq_class;
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "ZbasisComputation [inputMat] [output]\n";
      return -1;
    }
    // reading the matrix
    std::string FileI = argv[1];
    std::string FileO = argv[2];
    MyMatrix<T> TheMat = ReadMatrixFile<T>(FileI);
    // computing the kernel
    MyMatrix<T> TheBasis = GetZbasis(TheMat);
    //
    std::ofstream os(FileO);
    os << "return ";
    WriteMatrixGAP(os, TheBasis);
    os << ";\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
