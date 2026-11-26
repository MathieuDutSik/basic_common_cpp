// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm smoke test of the arithmetic functions of src_number on the boost
// types. ArithmeticFunctions.h and ResidueArithmetic.h are reached from
// COMB_Combinatorics.h, so most of polyhedral_common includes them: they
// must stay generic in the number type and not pull in GMP, which the Wasm
// builds do not have.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "ArithmeticFunctions.h"
#include "ResidueArithmetic.h"
#include "COMB_Combinatorics.h"
#include <iostream>

using cpp_int = boost::multiprecision::cpp_int;
using cpp_rational = boost::multiprecision::cpp_rational;

static int n_fail = 0;

#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << " : " << #cond    \
                << "\n";                                                       \
      n_fail++;                                                                \
    }                                                                          \
  } while (0)

static void test_arithmetic_functions() {
  std::optional<std::pair<cpp_int, size_t>> pp =
      PrimePowerDecomposition(cpp_int(343));
  CHECK(pp && pp->first == 7 && pp->second == 3);
  CHECK(!PrimePowerDecomposition(cpp_int(12)));
  CHECK((BinomialCoefficient<cpp_int, cpp_int>(cpp_int(10), cpp_int(3)) == 120));
  CHECK(NextPrimeAfter(13) == 17);
  CHECK(MultiplicativeOrder(2, 7) == 3);
  CHECK(IsFundamentalDiscriminant(-4));
  CHECK(!IsFundamentalDiscriminant(-16));
}

static void test_residue_arithmetic() {
  T_residue m7 = 7;
  CHECK(ResModOf(cpp_int(-1), m7) == 6);
  CHECK(ResInv(T_residue(3), m7) == 5);
  CHECK(ResPow(T_residue(2), uint64_t(10), T_residue(1000)) == 24);
  // 1/3 = 5 modulo 7.
  CHECK((ResModOfRational<cpp_rational, cpp_int>(cpp_rational(1, 3), m7) == 5));
}

int main() {
  test_arithmetic_functions();
  test_residue_arithmetic();
  if (n_fail > 0) {
    std::cerr << n_fail << " failure(s)\n";
    return 1;
  }
  std::cout << "Test_wasm_arithmetic: all checks passed\n";
  return 0;
}
