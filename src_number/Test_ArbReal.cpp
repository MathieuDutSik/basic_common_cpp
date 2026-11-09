// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#include "ArbReal.h"
#endif
// clang-format on

/*
  ArbReal: the enclosures of zeta(2) = pi^2 / 6 and Gamma(1/2)^2 = pi, and
  the round trips integer -> ball -> unique integer for the integer and
  rational types. Needs flint (ENABLE_FLINT_SUPPORT).
 */

#ifdef ENABLE_FLINT_SUPPORT
template <typename Tint, typename Tq> bool TestTypes(std::string const &name) {
  bool ok = true;
  for (std::string s : {"0", "-1", "123456789012345678901234567890",
                        "-98765432109876543210987654321"}) {
    Tint x = ParseScalar<Tint>(s);
    ArbReal b = arb_from_integer(x) + ArbReal(1) / ArbReal(3);
    b.add_error(ArbReal(1) / ArbReal(10));
    std::optional<Tint> y = arb_unique_integer<Tint>(b);
    if (y) {
      std::cerr << name << ": " << s << " + 1/3 +- 1/10 holds an integer\n";
      ok = false;
    }
    std::optional<Tint> z = arb_unique_integer<Tint>(arb_from_integer(x));
    if (!z || *z != x) {
      std::cerr << name << ": the round trip of " << s << " fails\n";
      ok = false;
    }
  }
  Tq q = ParseScalar<Tq>("-22/7");
  ArbReal r = arb_from_rational(q) * ArbReal(7);
  std::optional<Tint> m22 = arb_unique_integer<Tint>(r);
  if (!m22 || *m22 != ParseScalar<Tint>("-22")) {
    std::cerr << name << ": -22/7 * 7 gives " << r.str(20) << "\n";
    ok = false;
  }
  return ok;
}
#endif

int main() {
#ifdef ENABLE_FLINT_SUPPORT
  ArbPrecisionGuard g(200);
  bool ok = true;
  ArbReal pi = arb_pi();
  ArbReal d = arb_zeta(2) - pi * pi / ArbReal(6);
  ArbReal half = ArbReal(1) / ArbReal(2);
  ArbReal e = arb_gamma(half) * arb_gamma(half) - pi;
  if (!arb_contains_zero(d.v) || !arb_contains_zero(e.v) ||
      d.abs_ubound() > 1e-50 || e.abs_ubound() > 1e-50) {
    std::cerr << "zeta(2) or Gamma(1/2) wrong: " << d.str(10) << " "
              << e.str(10) << "\n";
    ok = false;
  }
  ok = TestTypes<mpz_class, mpq_class>("mpz_class") && ok;
  ok = TestTypes<boost::multiprecision::mpz_int,
                 boost::multiprecision::mpq_rational>("mpz_int") &&
       ok;
  ok = TestTypes<boost::multiprecision::cpp_int,
                 boost::multiprecision::cpp_rational>("cpp_int") &&
       ok;
  ok = TestTypes<fmpz_class, fmpq_class>("fmpz_class") && ok;
  std::cout << (ok ? "Test_ArbReal: all good" : "Test_ArbReal: FAILED") << "\n";
  return ok ? 0 : 1;
#else
  std::cout << "Test_ArbReal: needs flint (ENABLE_FLINT_SUPPORT), skipped\n";
  return 0;
#endif
}
