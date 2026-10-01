// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// Test of T_frexp, the mantissa-exponent form x = m 2^e, 1/2 <= |m| < 1, of
// every number type.
//
// The points checked are, for each type:
//   --- the exact dyadic cases, where m and e are known exactly: 0, 3 2^k,
//       -5 2^k, and for the fields 1 / (3 2^k);
//   --- values far beyond the range of a double, where a conversion to double
//       overflows but T_frexp must not: 3^700 and its inverse, to double
//       precision;
//   --- for QuadField, an element a + b sqrt(5) whose two terms cancel to
//       about 42 digits, which the conversion to double gets wrong and
//       T_frexp must get right through the conjugate;
//   --- for RealField, an element with coefficients of 400 digits.

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
// clang-format on
#include <cmath>
#include <string>
#include <vector>

static int n_error = 0;

static void check(bool test, std::string const &name) {
  if (test) {
    std::cerr << "PASS: " << name << "\n";
  } else {
    std::cerr << "FAIL: " << name << "\n";
    n_error++;
  }
}

// log2 |x| from the mantissa-exponent form.
template <typename T> double log2_abs(T const &x) {
  long e;
  double m = T_frexp(x, e);
  return std::log2(std::abs(m)) + static_cast<double>(e);
}

template <typename T> T power(T const &x, int k) {
  T ret(1);
  for (int i = 0; i < k; i++) {
    ret *= x;
  }
  return ret;
}

template <typename T>
bool is_exactly(T const &x, double const &m_expected, long const &e_expected) {
  long e;
  double m = T_frexp(x, e);
  return m == m_expected && e == e_expected;
}

template <typename T> bool is_normalized(T const &x) {
  long e;
  double m = T_frexp(x, e);
  return 0.5 <= std::abs(m) && std::abs(m) < 1.0;
}

// The checks that make sense for any type holding integers up to 2^k.
template <typename T> void test_integer(std::string const &name, int k) {
  T two(2);
  T pow2 = power(two, k);
  check(is_exactly<T>(T(0), 0.0, 0), name + ": 0");
  check(is_exactly<T>(T(1), 0.5, 1), name + ": 1");
  check(is_exactly<T>(T(3) * pow2, 0.75, k + 2), name + ": 3 2^k");
  check(is_exactly<T>(T(-5) * pow2, -0.625, k + 3), name + ": -5 2^k");
  check(is_normalized<T>(T(12345) * pow2), name + ": mantissa normalized");
}

// For the types holding arbitrarily large integers: 3^700, about 2^1109,
// beyond the range of a double.
template <typename T> void test_large_integer(std::string const &name) {
  test_integer<T>(name, 400);
  T x = power(T(3), 700);
  double expected = 700.0 * std::log2(3.0);
  check(std::abs(log2_abs<T>(x) - expected) < 1e-12, name + ": 3^700");
  check(std::abs(log2_abs<T>(T(-7) * x) - expected - std::log2(7.0)) < 1e-12,
        name + ": -7 3^700");
}

// For the fields: the same, and the inverses.
template <typename T> void test_large_rational(std::string const &name) {
  test_large_integer<T>(name);
  T pow2 = power(T(2), 400);
  check(is_exactly<T>(T(1) / (T(3) * pow2), 2.0 / 3.0, -401),
        name + ": 1 / (3 2^k)");
  T x = power(T(3), 700);
  double expected = -700.0 * std::log2(3.0);
  check(std::abs(log2_abs<T>(T(1) / x) - expected) < 1e-12, name + ": 3^-700");
  check(std::abs(log2_abs<T>(T(2) / T(3)) - std::log2(2.0 / 3.0)) < 1e-15,
        name + ": 2/3");
}

// QuadField<Trat, 5>: psi = (1 - sqrt(5)) / 2 has |psi| < 1, and psi^k has
// coefficients of the size of phi^k while its value is |psi|^k. At k = 200
// the two terms cancel to about 42 digits.
template <typename Trat> void test_quad_field(std::string const &name) {
  using T = QuadField<Trat, 5>;
  Trat half = Trat(1) / Trat(2);
  test_integer<T>(name, 400);
  check(is_exactly<T>(T(half, 0), 0.5, 0), name + ": 1/2");
  T phi(half, half);
  T psi(half, -half);
  int k = 200;
  double log2_phi = std::log2((1.0 + std::sqrt(5.0)) / 2.0);
  double log2_psi = std::log2((std::sqrt(5.0) - 1.0) / 2.0);
  check(std::abs(log2_abs<T>(power(phi, k)) - k * log2_phi) < 1e-11,
        name + ": phi^200, terms of the same sign");
  T psi_k = power(psi, k);
  check(std::abs(log2_abs<T>(psi_k) - k * log2_psi) < 1e-11,
        name + ": psi^200, terms cancelling");
  long e;
  double m = T_frexp(psi_k, e);
  check(m > 0, name + ": psi^200 is positive");
  T minus_psi_k = T(-1) * psi_k;
  check(std::abs(log2_abs<T>(minus_psi_k) - k * log2_psi) < 1e-11 &&
            T_frexp(minus_psi_k, e) < 0,
        name + ": -psi^200");
  // The conversion to double loses this value entirely, which is the point.
  double naive = UniversalScalarConversion<double, T>(psi_k);
  std::cerr << "  " << name << ": psi^200 through double is " << naive
            << " against " << std::pow(2.0, k * log2_psi) << "\n";
}

int const idx_field = 1;

// RealField over the cubic field of discriminant 49, the generator being
// 2 cos(2 pi / 7) of minimal polynomial X^3 + X^2 - 2X - 1.
void test_real_field() {
  using Tfield = RealField<idx_field>;
  using Tring = RealRing<idx_field>;
  using Tz = Tint_real_field;
  std::string eFile = "CI_tests/RealAlgebraicField/CubicFieldDisc_49";
  bool found = false;
  for (int level = 0; level <= 10; level++) {
    if (FILE_IsExistingFile(eFile)) {
      found = true;
      break;
    }
    eFile = "../" + eFile;
  }
  if (!found) {
    std::cerr << "Failed to find RealAlgebraicField test data\n";
    throw TerminalException{1};
  }
  HelperClassRealField<Trat_real_field> hcrf(eFile);
  insert_helper_real_algebraic_field(idx_field, hcrf);
  double theta = 2 * std::cos(2 * M_PI / 7);
  // 3^900 (1 + theta), with 3^900 about 2^1426.
  using Tq = Trat_real_field;
  Tq big = power(Tq(3), 900);
  Tfield x(std::vector<Tq>{big, big, Tq(0)});
  double expected = 900.0 * std::log2(3.0) + std::log2(1.0 + theta);
  check(std::abs(log2_abs<Tfield>(x) - expected) < 1e-12, "RealField: 3^900 (1+x)");
  check(std::abs(log2_abs<Tfield>(Tfield(1) / x) + expected) < 1e-12,
        "RealField: its inverse");
  Tz big_z = power(Tz(3), 900);
  Tring y(std::vector<Tz>{big_z, big_z, Tz(0)});
  check(std::abs(log2_abs<Tring>(y) - expected) < 1e-12, "RealRing: 3^900 (1+x)");
  check(is_exactly<Tfield>(Tfield(-6), -0.75, 3), "RealField: -6");
}

int main() {
  try {
    test_integer<int64_t>("int64_t", 40);
    test_integer<double>("double", 400);
    check(is_exactly<double>(0.1, 0.8, -3), "double: 0.1");
    test_integer<SafeInt64>("SafeInt64", 40);
    test_large_integer<mpz_class>("mpz_class");
    test_large_rational<mpq_class>("mpq_class");
    test_large_integer<boost::multiprecision::mpz_int>("mpz_int");
    test_large_rational<boost::multiprecision::mpq_rational>("mpq_rational");
    test_large_integer<boost::multiprecision::cpp_int>("cpp_int");
    test_large_rational<boost::multiprecision::cpp_rational>("cpp_rational");
    test_large_rational<Rational<mpz_class>>("Rational<mpz_class>");
#ifdef ENABLE_FLINT_SUPPORT
    test_large_integer<fmpz_class>("fmpz_class");
    test_large_rational<fmpq_class>("fmpq_class");
#endif
    test_quad_field<mpq_class>("QuadField<mpq_class, 5>");
    test_quad_field<boost::multiprecision::mpq_rational>(
        "QuadField<mpq_rational, 5>");
    test_real_field();
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
  if (n_error > 0) {
    std::cerr << "Test_Frexp: " << n_error << " failures\n";
    return 1;
  }
  std::cerr << "Test_Frexp: all checks pass\n";
  return 0;
}
