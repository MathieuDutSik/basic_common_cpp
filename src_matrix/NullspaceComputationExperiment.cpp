// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_Matrix.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  using T = mpq_class;
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "NullspaceComputation [inputMat] [output]\n";
      return -1;
    }
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> TheMat = ReadMatrix<T>(INmat);
    // computing the kernel
    Eigen::FullPivLU<MyMatrix<T>> lu(TheMat);
    MyMatrix<T> TheKer = lu.kernel();
    //
    //    MyMatrix<T> TheKer=NullspaceMat(TheMat);
    //
    std::ofstream os(argv[2]);
    os << "return ";
    WriteMatrixGAP(os, TheKer);
    os << ";\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
  runtime(time);
}
