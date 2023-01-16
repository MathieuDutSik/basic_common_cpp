// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "QuadField.h"

int main(int argc, char *argv[]) {
  using Trat = mpq_class;
  using T = QuadField<Trat, 5>;
  try {
    if (argc != 3) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "NullspaceComputationQ5 [inputMat] [output]\n");
      return -1;
    }
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> TheMat = ReadMatrix<T>(INmat);
    std::cerr << "TheMat=\n";
    WriteMatrix(std::cerr, TheMat);
    std::cerr << "\n";
    MyMatrix<T> TheKer = NullspaceMat(TheMat);
    std::cerr << "TheKer=\n";
    WriteMatrix(std::cerr, TheKer);
    std::cerr << "\n";
    //
    std::ofstream os(argv[2]);
    os << "return ";
    WriteMatrixGAP(os, TheKer);
    os << ";\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
