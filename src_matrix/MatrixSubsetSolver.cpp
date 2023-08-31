// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "MAT_Matrix_SubsetSolver.h"
// clang-format on


template<typename T>
void process(std::string const& FileI, std::string const& FileO) {
  using Tint = typename SubsetRankOneSolver<T>::Tint;
  MyMatrix<Tint> TheMat = ReadMatrixFile<Tint>(FileI);
  int n_row = TheMat.rows();
  // This is a test case
  SubsetRankOneSolver<T> solver(TheMat);
  for (int i_row=0; i_row<n_row; i_row++) {
    Face f(n_row);
    for (int j_row=0; j_row<n_row; j_row++) {
      if (i_row != j_row) {
        f[j_row] = 1;
      }
    }
    MyVector<Tint> V = solver.GetKernelVector(f);
    std::cerr << "i_row=" << i_row << " V=" << StringVector(V) << "\n";
  }
}


int main(int argc, char *argv[]) {
  try {
    if (argc != 4) {
      std::cerr << "This program is used as\n";
      std::cerr << "RowHermiteNormalForm [arith] [inputMat] [output]\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string FileO = argv[3];
    auto f=[&]() -> void {
      if (arith == "mpq_class") {
        using T = mpq_class;
        return process<T>(FileI, FileO);
      }
      if (arith == "cpp_rational") {
        using T = boost::multiprecision::cpp_rational;
        return process<T>(FileI, FileO);
      }
      if (arith == "safe_rational") {
        using T = Rational<SafeInt64>;
        return process<T>(FileI, FileO);
      }
      std::cerr << "Failed to find a matching type\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
}
