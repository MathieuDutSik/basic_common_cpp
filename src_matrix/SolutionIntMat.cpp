// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 2) {
      std::cerr << "This program is used as\n";
      std::cerr << "SolutionIntMat [inputMat] [output]\n";
      std::cerr << "  or\n";
      std::cerr << "SolutionIntMat [inputMat]\n";
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
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
