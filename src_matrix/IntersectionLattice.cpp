// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "rational.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "QuadField.h"
#include "MAT_MatrixInt.h"
// clang-format off

template<typename T>
void process() {
  int n = 4;
  auto get_fullrank = [&]() -> MyMatrix<T> {
    MyMatrix<T> M(n, n);
    while (true) {
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
          M(i, j) = random() % 7;
      if (RankMat(M) == n)
        return M;
    }
  };
  for (int iter = 0; iter < 10; iter++) {
    std::cerr << "iter=" << iter << "\n";
    MyMatrix<T> M1 = get_fullrank();
    MyMatrix<T> M2 = get_fullrank();
    MyMatrix<T> M1_i_M2 = IntersectionLattice(M1, M2);
    for (int i = 0; i < n; i++) {
      MyVector<T> v = GetMatrixRow(M1_i_M2, i);
      if (!SolutionIntMat(M1, v) || !SolutionIntMat(M2, v)) {
        std::cerr << "   M1=\n";
        WriteMatrix(std::cerr, M1);
        std::cerr << "   M2=\n";
        WriteMatrix(std::cerr, M2);
        std::cerr << "v=";
        WriteVector(std::cerr, v);
        std::cerr << "Could not resolve v in M1 or M2\n";
        throw TerminalException{1};
      }
    }
  }
}



int main(int argc, char *argv[]) {
  //  using T=mpz_class;
  //  using T=mpq_class;
  using T = Rational<long long>;
  try {
    if (argc != 2) {
      std::cerr << "IntersectionLattice [arithmetic]\n";
      throw TerminalException{1};
    }
    std::string arith = argv[1];
    auto f=[&]() -> void {
      if (arith == "mpz_class")
        return process<mpz_class>();
      if (arith == "mpq_class")
        return process<mpq_class>();
      if (arith == "Rational<long long>")
        return process<Rational<long long>>();
      //      if (arith == "safe_integer")
      //        return process<SafeInt64>();
      //      if (arith == "safe_rational")
      //        return process<Rational<SafeInt64>>();
      std::cerr << "Failed to have a matching entry\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
    std::cerr << "Erroneous termination of the program\n";
  }
}
