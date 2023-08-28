// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MatrixLinbox.h"
// clang-format off

int main(int argc, char *argv[]) {
  using T = mpq_class;
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "NullspaceComputationLinbox [inputMat] [output]\n";
      return -1;
    }
    // reading the matrix
    std::string eFileI = argv[1];
    std::string eFileO = argv[2];
    MyMatrix<T> TheMat = ReadMatrixFile<T>(eFileI);
    MyMatrix<T> TheKer = NullspaceTrMat_linbox(TheMat);
    //
    std::ofstream os(eFileO);
    os << "return ";
    WriteMatrixGAP(os, TheKer);
    os << ";\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
