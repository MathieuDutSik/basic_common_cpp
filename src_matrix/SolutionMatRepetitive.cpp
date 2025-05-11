// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_Matrix.h"
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
    using T = mpq_class;
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> TheMat = ReadMatrix<T>(INmat);
    MyVector<T> TheVec = ReadVector<T>(INmat);
    SolutionMatRepetitive<T> smr(TheMat);
    std::optional<MyVector<T>> result = smr.GetSolution(TheVec);
    // computing one solution (maybe)
    auto Prt = [&](std::ostream &os) -> void {
      if (result) {
        os << "return " << StringVectorGAP(*result) << ";\n";
      } else {
        os << "return fail;\n";
      }
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
