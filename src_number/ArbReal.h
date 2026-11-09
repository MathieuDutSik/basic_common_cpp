// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_ARBREAL_H_
#define SRC_NUMBER_ARBREAL_H_

// clang-format off
#include "Basic_functions.h"
#include <flint/arb.h>
#include <flint/arb_hypgeom.h>
#include <flint/fmpz.h>
#include <optional>
#include <sstream>
#include <string>
// clang-format on

/*
  ArbReal: a real ball of Arb (flint), with value semantics. The operations
  use the working precision ArbReal::prec(), set with ArbPrecisionGuard;
  every result is a ball containing the exact result of the operation on
  the points of the input balls, so a computation made of these operations
  gives rigorous enclosures.

  The conversions from and to the exact types go through the decimal
  strings, for any integer type (int64_t, mpz_class, fmpz_class, boost
  mpz_int, cpp_int) and any rational type printed as "a" or "a/b".
 */
class ArbReal {
public:
  arb_t v;
  static slong &prec() {
    thread_local slong p = 64;
    return p;
  }
  ArbReal() { arb_init(v); }
  ArbReal(slong x) {
    arb_init(v);
    arb_set_si(v, x);
  }
  ArbReal(ArbReal const &x) {
    arb_init(v);
    arb_set(v, x.v);
  }
  ArbReal(ArbReal &&x) noexcept {
    arb_init(v);
    arb_swap(v, x.v);
  }
  ArbReal &operator=(ArbReal const &x) {
    arb_set(v, x.v);
    return *this;
  }
  ArbReal &operator=(ArbReal &&x) noexcept {
    arb_swap(v, x.v);
    return *this;
  }
  ~ArbReal() { arb_clear(v); }
  ArbReal &operator+=(ArbReal const &x) {
    arb_add(v, v, x.v, prec());
    return *this;
  }
  ArbReal &operator-=(ArbReal const &x) {
    arb_sub(v, v, x.v, prec());
    return *this;
  }
  ArbReal &operator*=(ArbReal const &x) {
    arb_mul(v, v, x.v, prec());
    return *this;
  }
  ArbReal &operator/=(ArbReal const &x) {
    arb_div(v, v, x.v, prec());
    return *this;
  }
  // the radius enlarged by e >= 0 (an upper bound)
  void add_error(ArbReal const &e) { arb_add_error(v, e.v); }
  // an upper bound of |x| as a double, rounded up
  double abs_ubound() const {
    arf_t u;
    arf_init(u);
    arb_get_abs_ubound_arf(u, v, 53);
    double d = arf_get_d(u, ARF_RND_UP);
    arf_clear(u);
    return d;
  }
  double mid_d() const { return arf_get_d(arb_midref(v), ARF_RND_NEAR); }
  bool is_positive() const { return arb_is_positive(v); }
  std::string str(slong digits) const {
    char *s = arb_get_str(v, digits, 0);
    std::string res(s);
    flint_free(s);
    return res;
  }
};

class ArbPrecisionGuard {
  slong old;

public:
  explicit ArbPrecisionGuard(slong p) : old(ArbReal::prec()) {
    ArbReal::prec() = p;
  }
  ~ArbPrecisionGuard() { ArbReal::prec() = old; }
};

inline ArbReal operator+(ArbReal x, ArbReal const &y) { return x += y; }
inline ArbReal operator-(ArbReal x, ArbReal const &y) { return x -= y; }
inline ArbReal operator*(ArbReal x, ArbReal const &y) { return x *= y; }
inline ArbReal operator/(ArbReal x, ArbReal const &y) { return x /= y; }
inline ArbReal operator-(ArbReal const &x) {
  ArbReal r;
  arb_neg(r.v, x.v);
  return r;
}

inline ArbReal arb_exp(ArbReal const &x) {
  ArbReal r;
  arb_exp(r.v, x.v, ArbReal::prec());
  return r;
}
inline ArbReal arb_log(ArbReal const &x) {
  ArbReal r;
  arb_log(r.v, x.v, ArbReal::prec());
  return r;
}
inline ArbReal arb_abs(ArbReal const &x) {
  ArbReal r;
  arb_abs(r.v, x.v);
  return r;
}
inline ArbReal arb_pow(ArbReal const &x, ArbReal const &y) {
  ArbReal r;
  arb_pow(r.v, x.v, y.v, ArbReal::prec());
  return r;
}
inline ArbReal arb_pow_si(ArbReal const &x, slong e) {
  ArbReal r;
  if (e >= 0) {
    arb_pow_ui(r.v, x.v, static_cast<ulong>(e), ArbReal::prec());
  } else {
    arb_pow_ui(r.v, x.v, static_cast<ulong>(-e), ArbReal::prec());
    arb_inv(r.v, r.v, ArbReal::prec());
  }
  return r;
}
inline ArbReal arb_cos(ArbReal const &x) {
  ArbReal r;
  arb_cos(r.v, x.v, ArbReal::prec());
  return r;
}
inline ArbReal arb_gamma(ArbReal const &x) {
  ArbReal r;
  arb_hypgeom_gamma(r.v, x.v, ArbReal::prec());
  return r;
}
inline ArbReal arb_pi() {
  ArbReal r;
  arb_const_pi(r.v, ArbReal::prec());
  return r;
}
inline ArbReal arb_euler() {
  ArbReal r;
  arb_const_euler(r.v, ArbReal::prec());
  return r;
}
inline ArbReal arb_zeta(ulong k) {
  ArbReal r;
  arb_zeta_ui(r.v, k, ArbReal::prec());
  return r;
}
// the ball [x - e, x + e] for a double x and e >= 0
inline ArbReal arb_from_double(double x) {
  ArbReal r;
  arb_set_d(r.v, x);
  return r;
}

template <typename Tint> ArbReal arb_from_integer(Tint const &x) {
  std::stringstream ss;
  ss << x;
  ArbReal r;
  arb_set_str(r.v, ss.str().c_str(), ArbReal::prec());
  return r;
}

template <typename T> ArbReal arb_from_rational(T const &x) {
  std::stringstream ss;
  ss << x;
  std::string s = ss.str();
  size_t pos = s.find('/');
  ArbReal r;
  if (pos == std::string::npos) {
    arb_set_str(r.v, s.c_str(), ArbReal::prec());
    return r;
  }
  ArbReal den;
  arb_set_str(r.v, s.substr(0, pos).c_str(), ArbReal::prec());
  arb_set_str(den.v, s.substr(pos + 1).c_str(), ArbReal::prec());
  return r / den;
}

// The integer in the ball, when there is exactly one.
template <typename Tint>
std::optional<Tint> arb_unique_integer(ArbReal const &x) {
  fmpz_t z;
  fmpz_init(z);
  int ok = arb_get_unique_fmpz(z, x.v);
  std::optional<Tint> res;
  if (ok) {
    char *s = fmpz_get_str(nullptr, 10, z);
    res = ParseScalar<Tint>(std::string(s));
    flint_free(s);
  }
  fmpz_clear(z);
  return res;
}

// clang-format off
#endif  // SRC_NUMBER_ARBREAL_H_
// clang-format on
