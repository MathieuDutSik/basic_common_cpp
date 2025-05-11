// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "MAT_Matrix_SubsetSolver.h"
// clang-format on

template <typename T>
inline
    typename std::enable_if<!has_reduction_subset_solver<T>::value, void>::type
    process(std::string const &FileI, [[maybe_unused]] std::string const &FileO) {
  MyMatrix<T> TheMat = ReadMatrixFile<T>(FileI);
  int n_row = TheMat.rows();
  int n_col = TheMat.cols();
  // This is a test case
  SubsetRankOneSolver_Field<T> solver(TheMat);
  for (int i_row = 0; i_row < n_row; i_row++) {
    Face f(n_row);
    for (int j_row = 0; j_row < n_row; j_row++) {
      if (i_row != j_row) {
        f[j_row] = 1;
      }
    }
    MyVector<T> V = solver.GetKernelVector(f);
    for (int j_row = 0; j_row < n_row; j_row++) {
      if (j_row != i_row) {
        T scal(0);
        for (int i_col = 0; i_col < n_col; i_col++) {
          scal += V(i_col) * TheMat(j_row, i_col);
        }
        if (scal != 0) {
          std::cerr << "j_row=" << j_row << " V=" << StringVector(V) << "\n";
          std::cerr << "Wrong scalar product\n";
          throw TerminalException{1};
        }
      }
    }
    std::cerr << "SCH A : i_row=" << i_row << " V=" << StringVector(V) << "\n";
  }
}

template <typename T>
inline
    typename std::enable_if<has_reduction_subset_solver<T>::value, void>::type
    process(std::string const &FileI, [[maybe_unused]] std::string const &FileO) {
  using Tint = typename SubsetRankOneSolver_Acceleration<T>::Tint;
  MyMatrix<Tint> TheMat = ReadMatrixFile<Tint>(FileI);
  int n_row = TheMat.rows();
  int n_col = TheMat.cols();
  // This is a test case
  SubsetRankOneSolver_Acceleration<T> solver(TheMat);
  for (int i_row = 0; i_row < n_row; i_row++) {
    Face f(n_row);
    for (int j_row = 0; j_row < n_row; j_row++) {
      if (i_row != j_row) {
        f[j_row] = 1;
      }
    }
    MyVector<Tint> V = solver.GetKernelVector(f);
    for (int j_row = 0; j_row < n_row; j_row++) {
      if (j_row != i_row) {
        Tint scal(0);
        for (int i_col = 0; i_col < n_col; i_col++) {
          scal += V(i_col) * TheMat(j_row, i_col);
        }
        if (scal != 0) {
          std::cerr << "j_row=" << j_row << " V=" << StringVector(V) << "\n";
          std::cerr << "Wrong scalar product\n";
          throw TerminalException{1};
        }
      }
    }
    std::cerr << "SCH B : i_row=" << i_row << " V=" << StringVector(V) << "\n";
  }
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4) {
      std::cerr << "This program is used as\n";
      std::cerr << "RowHermiteNormalForm [arith] [inputMat] [output]\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string FileO = argv[3];
    auto f = [&]() -> void {
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
      std::cerr << "Allowed types are mpq_class, cpp_rational, safe_rational\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
