// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_MPREALCONVERSION_H_
#define SRC_NUMBER_MPREALCONVERSION_H_

// mpreal 3.6.8 calls mpfr_get_prec on a temporary, which the macro version
// of MPFR 4 refuses
#ifndef MPFR_USE_NO_MACRO
#define MPFR_USE_NO_MACRO
#endif

// clang-format off
#include "Basic_functions.h"
#include "mpreal.h"
#include <sstream>
#include <string>
// clang-format on

/*
  Conversions between the exact number types and mpfr::mpreal, for any
  integer type (int64_t, mpz_class, fmpz_class, boost mpz_int, cpp_int) and
  any rational type printed as "a" or "a/b" (mpq_class, fmpq_class, boost
  mpq_rational, cpp_rational). They go through the decimal strings, so that
  no call specific to a library is needed.
 */

template <typename Tint>
mpfr::mpreal mpreal_from_integer(Tint const &x, mpfr_prec_t prec) {
  std::stringstream ss;
  ss << x;
  return mpfr::mpreal(ss.str(), prec, 10, MPFR_RNDN);
}

template <typename T>
mpfr::mpreal mpreal_from_rational(T const &x, mpfr_prec_t prec) {
  std::stringstream ss;
  ss << x;
  std::string s = ss.str();
  size_t pos = s.find('/');
  if (pos == std::string::npos) {
    return mpfr::mpreal(s, prec, 10, MPFR_RNDN);
  }
  mpfr::mpreal num(s.substr(0, pos), prec, 10, MPFR_RNDN);
  mpfr::mpreal den(s.substr(pos + 1), prec, 10, MPFR_RNDN);
  return num / den;
}

// The integer nearest to x, the halves rounded away from zero.
template <typename Tint> Tint mpreal_nearest_integer(mpfr::mpreal const &x) {
  mpfr::mpreal y = mpfr::round(x);
  char *s = nullptr;
  mpfr_asprintf(&s, "%.0Rf", y.mpfr_srcptr());
  std::string str(s);
  mpfr_free_str(s);
  return ParseScalar<Tint>(str);
}

// clang-format off
#endif  // SRC_NUMBER_MPREALCONVERSION_H_
// clang-format on
