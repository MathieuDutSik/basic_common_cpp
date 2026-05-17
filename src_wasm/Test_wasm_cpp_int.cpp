// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// First Wasm smoke test: exercises boost::multiprecision::cpp_int and
// cpp_rational arithmetic. No GMP, no threads, no fork/exec — all the
// pieces split out in Basic_external_program.h / Basic_threading.h are
// kept out of the include graph.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "TypeConversion.h"
#include "TypeConversionFinal.h"
#include "factorizations.h"
#include "fractions.h"
#include "quadratic_residue.h"
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

static void test_cpp_int_arithmetic() {
  cpp_int a("123456789012345678901234567890");
  cpp_int b("987654321098765432109876543210");
  cpp_int sum = a + b;
  cpp_int prod = a * b;
  CHECK(sum == cpp_int("1111111110111111111011111111100"));
  CHECK(prod / b == a);
  CHECK(prod % b == 0);
  CHECK(a - a == 0);
  CHECK((a * 0) == 0);
}

static void test_cpp_rational_arithmetic() {
  cpp_rational a(1, 3);
  cpp_rational b(1, 6);
  cpp_rational sum = a + b;
  CHECK(sum == cpp_rational(1, 2));
  cpp_rational prod = a * b;
  CHECK(prod == cpp_rational(1, 18));
  cpp_rational quot = a / b;
  CHECK(quot == cpp_rational(2));
}

static void test_factorization_cpp_int() {
  cpp_int n("12345678901234567890");
  // 12345678901234567890 = 2 * 3 * 3 * 5 * 101 * 3541 * 3607 * 27961 * 3803
  std::map<cpp_int, size_t> factors = FactorsIntMap(n);
  cpp_int recomposed = 1;
  for (auto const &kv : factors) {
    for (size_t i = 0; i < kv.second; i++) {
      recomposed *= kv.first;
    }
  }
  CHECK(recomposed == n);
}

static void test_quadratic_residue_cpp_int() {
  // 23 is prime, so exactly (23-1)/2 + 1 = 12 residues in [0, 22].
  cpp_int p(23);
  int n_quad = 0;
  for (int u = 0; u < 23; u++) {
    cpp_int u_T(u);
    if (is_quadratic_residue(u_T, p)) {
      n_quad++;
    }
  }
  CHECK(n_quad == 12);
}

static void test_floor_ceil_cpp_rational() {
  cpp_rational x = cpp_rational(7, 3);
  cpp_int floor_x = UniversalFloorScalarInteger<cpp_int, cpp_rational>(x);
  cpp_int ceil_x = UniversalCeilScalarInteger<cpp_int, cpp_rational>(x);
  CHECK(floor_x == 2);
  CHECK(ceil_x == 3);
  cpp_rational y = cpp_rational(-7, 3);
  cpp_int floor_y = UniversalFloorScalarInteger<cpp_int, cpp_rational>(y);
  cpp_int ceil_y = UniversalCeilScalarInteger<cpp_int, cpp_rational>(y);
  CHECK(floor_y == -3);
  CHECK(ceil_y == -2);
}

static void test_continued_fraction_cpp_rational() {
  // 22/7 has known continued fraction [3; 7].
  cpp_rational x = cpp_rational(22, 7);
  std::vector<cpp_rational> seq =
      get_sequence_continuous_fraction_approximant(x);
  CHECK(!seq.empty());
  // Last entry of the sequence equals x.
  CHECK(seq.back() == x);
}

int main() {
  test_cpp_int_arithmetic();
  test_cpp_rational_arithmetic();
  test_factorization_cpp_int();
  test_quadratic_residue_cpp_int();
  test_floor_ceil_cpp_rational();
  test_continued_fraction_cpp_rational();
  if (n_fail == 0) {
    std::cerr << "Test_wasm_cpp_int: all checks passed.\n";
    return 0;
  }
  std::cerr << "Test_wasm_cpp_int: " << n_fail << " check(s) failed.\n";
  return 1;
}
