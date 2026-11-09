// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "MprealConversion.h"
// clang-format on

/*
  The round trip integer -> mpreal -> nearest integer, and the rationals,
  for the integer and rational types.
 */

template <typename Tint, typename Tq> bool TestTypes(std::string const &name) {
  bool ok = true;
  std::vector<std::string> ints{"0", "1", "-1",
                                "123456789012345678901234567890",
                                "-98765432109876543210987654321"};
  for (auto const &s : ints) {
    Tint x = ParseScalar<Tint>(s);
    mpfr::mpreal y = mpreal_from_integer(x, 256);
    Tint z = mpreal_nearest_integer<Tint>(y + mpfr::mpreal("0.3", 256));
    Tint w = mpreal_nearest_integer<Tint>(y - mpfr::mpreal("0.3", 256));
    if (z != x || w != x) {
      std::cerr << name << ": round trip of " << s << " gives " << z << " " << w
                << "\n";
      ok = false;
    }
  }
  Tq q = ParseScalar<Tq>("-22/7");
  mpfr::mpreal r = mpreal_from_rational(q, 128);
  mpfr::mpreal ref = mpfr::mpreal(-22, 128) / 7;
  if (abs(r - ref) > mpfr::mpreal("1e-35", 128)) {
    std::cerr << name << ": -22/7 gives " << r << "\n";
    ok = false;
  }
  Tq q2 = ParseScalar<Tq>("5");
  if (mpreal_from_rational(q2, 64) != 5) {
    std::cerr << name << ": 5 gives " << mpreal_from_rational(q2, 64) << "\n";
    ok = false;
  }
  if (mpreal_nearest_integer<Tint>(mpfr::mpreal("-2.5", 64)) !=
      ParseScalar<Tint>("-3")) {
    std::cerr << name << ": -2.5 is not rounded to -3\n";
    ok = false;
  }
  return ok;
}

int main() {
  bool ok = true;
  ok = TestTypes<mpz_class, mpq_class>("mpz_class") && ok;
  ok = TestTypes<boost::multiprecision::mpz_int,
                 boost::multiprecision::mpq_rational>("mpz_int") &&
       ok;
  ok = TestTypes<boost::multiprecision::cpp_int,
                 boost::multiprecision::cpp_rational>("cpp_int") &&
       ok;
#ifdef ENABLE_FLINT_SUPPORT
  ok = TestTypes<fmpz_class, fmpq_class>("fmpz_class") && ok;
#endif
  std::cout << (ok ? "Test_MprealConversion: all good"
                   : "Test_MprealConversion: FAILED")
            << "\n";
  return ok ? 0 : 1;
}
