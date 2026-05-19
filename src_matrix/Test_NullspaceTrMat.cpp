// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "MAT_Matrix.h"
// clang-format on

template <typename T> void test_type() {
  MyMatrix<T> M = ZeroMatrix<T>(3, 4);
  MyMatrix<T> NSP = NullspaceTrMat(M);
  std::cerr << "rows=" << NSP.rows() << " cols=" << NSP.cols() << "\n";
}

int main() {
  try {
    test_type<mpq_class>();
    test_type<mpz_class>();
    test_type<boost::multiprecision::cpp_int>();
    test_type<boost::multiprecision::cpp_rational>();
    test_type<boost::multiprecision::mpz_int>();
    test_type<boost::multiprecision::mpq_rational>();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
  return 0;
}
