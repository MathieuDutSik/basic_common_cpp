// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Eigen-compatibility compile test.
//
// Verifies that MyMatrix<T> matrix arithmetic compiles against the linked
// Eigen for every scalar type the project supports:
//
//   built-in:   int / SafeInt64 / Rational<SafeInt64>
//   GMP-C++:    mpz_class / mpq_class
//   boost mp:   cpp_int / cpp_rational / mpz_int / mpq_rational
//
// The boost-multiprecision <-> Eigen adapter is opted into here, which in
// turn pulls in src_matrix/EigenBoostNumTraits.h (the Literal=int override
// that lets Eigen >= 3.5 work with boost cpp_int / cpp_rational /
// mpz_int / mpq_rational without comparing against literal 0.0).
//
// The matrix multiplications below all trigger Eigen's is_exactly_zero
// path internally; if any of the integer/rational backends regresses, this
// test will fail to compile.
//
// This is a compile + smoke test only — each multiplication's result is
// trivially the identity matrix, which the test asserts.

// clang-format off
#define INCLUDE_NUMBER_THEORY_BOOST_CPP_INT
#define INCLUDE_NUMBER_THEORY_BOOST_GMP_INT

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "MAT_Matrix.h"
// clang-format on

#include <iostream>

template <typename T> bool check_identity_product(std::string const &name) {
  int n = 4;
  MyMatrix<T> I = IdentityMat<T>(n);
  MyMatrix<T> P = I * I;
  if (!TestEqualityMatrix(P, I)) {
    std::cerr << "FAIL: I*I != I for " << name << "\n";
    return false;
  }
  return true;
}

int main() {
  try {
    bool ok = true;
    ok &= check_identity_product<mpz_class>("mpz_class");
    ok &= check_identity_product<mpq_class>("mpq_class");
    ok &= check_identity_product<boost::multiprecision::cpp_int>("cpp_int");
    ok &= check_identity_product<boost::multiprecision::cpp_rational>(
        "cpp_rational");
    ok &= check_identity_product<boost::multiprecision::mpz_int>("mpz_int");
    ok &= check_identity_product<boost::multiprecision::mpq_rational>(
        "mpq_rational");
    ok &= check_identity_product<SafeInt64>("SafeInt64");
    ok &= check_identity_product<Rational<SafeInt64>>("Rational<SafeInt64>");
    if (!ok) {
      std::cerr << "Test_eigen_all_types: at least one type failed\n";
      return 1;
    }
    std::cerr << "Test_eigen_all_types: all 8 types compiled and ran\n";
    return 0;
  } catch (TerminalException const &e) {
    std::cerr << "Test_eigen_all_types: TerminalException\n";
    return e.eVal;
  }
}
