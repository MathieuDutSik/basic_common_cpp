// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format on

int main(int argc, char *argv[]) {
  try {
    if (argc != 3 && argc != 2) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "SolutionIntMat [inputMat] [output]\n");
      fprintf(stderr, "or\n");
      fprintf(stderr, "SolutionIntMat [inputMat]\n");
      return -1;
    }
    using Tint = mpz_class;
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<Tint> TheMat = ReadMatrix<Tint>(INmat);
    MyVector<Tint> TheVec = ReadVector<Tint>(INmat);
    std::optional<MyVector<Tint>> result = SolutionIntMat(TheMat, TheVec);
    // computing one solution (maybe)
    auto Prt = [&](std::ostream &os) -> void {
      os << "return " + ResultSolutionIntMat_to_GAP(result) + ";\n";
    };
    //
    if (argc == 2) {
      Prt(std::cerr);
    } else {
      std::ofstream os(argv[2]);
      Prt(os);
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
