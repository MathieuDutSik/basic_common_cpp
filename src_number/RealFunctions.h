// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_REALFUNCTIONS_H_
#define SRC_NUMBER_REALFUNCTIONS_H_

// clang-format off
#include "DoubleDouble.h"
#include "MprealConversion.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "ArbReal.h"
#endif
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <type_traits>
// clang-format on

/*
  One interface to the real types double, DoubleDouble, mpfr::mpreal and
  ArbReal (the balls of Arb, part of FLINT, hence only with
  ENABLE_FLINT_SUPPORT), so that a numerical computation can be written
  once as a template on Tfloat:
  - the constants real_pi, real_euler, real_zeta(k) and real_infinity, and
    real_from_double, real_from_int, real_from_integer (any integer type);
  - real_exp, real_log, real_cos, real_abs, real_gamma, real_pow_int, and
    real_to_double (the midpoint of a ball);
  - RealPrecisionGuard<Tfloat>(bits): the working precision of mpreal and
    ArbReal while it lives, nothing for double and DoubleDouble;
  - RealIsBall<Tfloat>: whether the type encloses its errors (ArbReal), and
    then real_add_error(x, b) widens the ball x by b, real_upper(x) is an
    exact upper bound of |x| and real_certainly_le(x, y) holds when every
    point of x is at most every point of y; for the other types they are
    the identity, |x| and x <= y;
  - real_convert<Tto>(x): x in another type, at the working precision of
    Tto.
  The transcendental functions of DoubleDouble beyond exp and log (cos,
  gamma, zeta and the constants) go through mpreal at 128 bits.
 */

template <typename T> struct RealIsBallTrait {
  static constexpr bool value = false;
};
#ifdef ENABLE_FLINT_SUPPORT
template <> struct RealIsBallTrait<ArbReal> {
  static constexpr bool value = true;
};
#endif
template <typename T> constexpr bool RealIsBall = RealIsBallTrait<T>::value;

// The precision guards
template <typename T> class RealPrecisionGuard {
public:
  explicit RealPrecisionGuard(long) {}
};

template <> class RealPrecisionGuard<mpfr::mpreal> {
  mpfr_prec_t old;

public:
  explicit RealPrecisionGuard(long bits) : old(mpfr::mpreal::get_default_prec()) {
    mpfr::mpreal::set_default_prec(static_cast<mpfr_prec_t>(bits));
  }
  ~RealPrecisionGuard() { mpfr::mpreal::set_default_prec(old); }
  RealPrecisionGuard(RealPrecisionGuard const &) = delete;
  RealPrecisionGuard &operator=(RealPrecisionGuard const &) = delete;
};

#ifdef ENABLE_FLINT_SUPPORT
template <> class RealPrecisionGuard<ArbReal> {
  ArbPrecisionGuard g;

public:
  explicit RealPrecisionGuard(long bits) : g(static_cast<slong>(bits)) {}
};
#endif

// double
inline double real_exp(double x) { return std::exp(x); }
inline double real_log(double x) { return std::log(x); }
inline double real_cos(double x) { return std::cos(x); }
inline double real_abs(double x) { return std::abs(x); }
inline double real_gamma(double x) { return std::tgamma(x); }
inline double real_to_double(double x) { return x; }
inline double real_upper(double x) { return std::abs(x); }
inline void real_add_error(double &, double const &) {}
inline bool real_certainly_le(double x, double y) { return x <= y; }

// mpreal
inline mpfr::mpreal real_exp(mpfr::mpreal const &x) { return mpfr::exp(x); }
inline mpfr::mpreal real_log(mpfr::mpreal const &x) { return mpfr::log(x); }
inline mpfr::mpreal real_cos(mpfr::mpreal const &x) { return mpfr::cos(x); }
inline mpfr::mpreal real_abs(mpfr::mpreal const &x) { return mpfr::abs(x); }
inline mpfr::mpreal real_gamma(mpfr::mpreal const &x) {
  return mpfr::gamma(x);
}
inline double real_to_double(mpfr::mpreal const &x) { return x.toDouble(); }
inline mpfr::mpreal real_upper(mpfr::mpreal const &x) { return mpfr::abs(x); }
inline void real_add_error(mpfr::mpreal &, mpfr::mpreal const &) {}
inline bool real_certainly_le(mpfr::mpreal const &x, mpfr::mpreal const &y) {
  return x <= y;
}

// DoubleDouble
inline DoubleDouble DoubleDoubleFromMpreal(mpfr::mpreal const &x) {
  double hi = x.toDouble();
  double lo = (x - mpfr::mpreal(hi, x.get_prec())).toDouble();
  return DoubleDoubleQuickTwoSum(hi, lo);
}
inline mpfr::mpreal DoubleDoubleToMpreal(DoubleDouble const &x,
                                         mpfr_prec_t prec) {
  return mpfr::mpreal(x.hi, prec) + mpfr::mpreal(x.lo, prec);
}
inline DoubleDouble real_exp(DoubleDouble const &x) {
  return DoubleDoubleExp(x);
}
inline DoubleDouble real_log(DoubleDouble const &x) {
  return DoubleDoubleLog(x);
}
inline DoubleDouble real_cos(DoubleDouble const &x) {
  return DoubleDoubleFromMpreal(mpfr::cos(DoubleDoubleToMpreal(x, 128)));
}
inline DoubleDouble real_abs(DoubleDouble const &x) {
  return (x.hi < 0) ? -x : x;
}
inline DoubleDouble real_gamma(DoubleDouble const &x) {
  return DoubleDoubleFromMpreal(mpfr::gamma(DoubleDoubleToMpreal(x, 128)));
}
inline double real_to_double(DoubleDouble const &x) { return x.hi + x.lo; }
inline DoubleDouble real_upper(DoubleDouble const &x) { return real_abs(x); }
inline void real_add_error(DoubleDouble &, DoubleDouble const &) {}
inline bool real_certainly_le(DoubleDouble const &x, DoubleDouble const &y) {
  return x <= y;
}

#ifdef ENABLE_FLINT_SUPPORT
// ArbReal
inline ArbReal real_exp(ArbReal const &x) { return arb_exp(x); }
inline ArbReal real_log(ArbReal const &x) { return arb_log(x); }
inline ArbReal real_cos(ArbReal const &x) { return arb_cos(x); }
inline ArbReal real_abs(ArbReal const &x) { return arb_abs(x); }
inline ArbReal real_gamma(ArbReal const &x) { return arb_gamma(x); }
inline double real_to_double(ArbReal const &x) { return x.mid_d(); }
// an exact upper bound of |x|: a ball of radius 0
inline ArbReal real_upper(ArbReal const &x) {
  arf_t u;
  arf_init(u);
  arb_get_abs_ubound_arf(u, x.v, ArbReal::prec());
  ArbReal r;
  arb_set_arf(r.v, u);
  arf_clear(u);
  return r;
}
inline void real_add_error(ArbReal &x, ArbReal const &b) { x.add_error(b); }
inline bool real_certainly_le(ArbReal const &x, ArbReal const &y) {
  return arb_le(x.v, y.v);
}
#endif

// The constants, at the working precision
template <typename T> T real_from_double(double x);
template <> inline double real_from_double<double>(double x) { return x; }
template <> inline DoubleDouble real_from_double<DoubleDouble>(double x) {
  return DoubleDouble(x);
}
template <>
inline mpfr::mpreal real_from_double<mpfr::mpreal>(double x) {
  return mpfr::mpreal(x);
}
#ifdef ENABLE_FLINT_SUPPORT
template <> inline ArbReal real_from_double<ArbReal>(double x) {
  return arb_from_double(x);
}
#endif

template <typename T> T real_from_int(int64_t x) {
  if constexpr (std::is_same_v<T, mpfr::mpreal>) {
    return mpfr::mpreal(static_cast<long>(x));
  } else if constexpr (std::is_same_v<T, DoubleDouble>) {
    return DoubleDouble(x);
  } else if constexpr (std::is_same_v<T, double>) {
    return static_cast<double>(x);
  } else {
    return T(static_cast<long>(x));
  }
}

// An integer of any type Tint: exactly in a ball, rounded otherwise.
template <typename T, typename Tint> T real_from_integer(Tint const &x) {
#ifdef ENABLE_FLINT_SUPPORT
  if constexpr (RealIsBall<T>) {
    return arb_from_integer(x);
  } else
#endif
  if constexpr (std::is_same_v<T, mpfr::mpreal>) {
    return mpreal_from_integer(x, mpfr::mpreal::get_default_prec());
  } else if constexpr (std::is_same_v<T, DoubleDouble>) {
    return DoubleDoubleFromMpreal(mpreal_from_integer(x, 128));
  } else {
    return mpreal_from_integer(x, 64).toDouble();
  }
}

template <typename T> T real_pi() {
  if constexpr (std::is_same_v<T, double>) {
    return M_PI;
  } else if constexpr (std::is_same_v<T, DoubleDouble>) {
    return DoubleDouble(3.141592653589793, 1.2246467991473532e-16);
  } else if constexpr (std::is_same_v<T, mpfr::mpreal>) {
    return mpfr::const_pi();
  }
#ifdef ENABLE_FLINT_SUPPORT
  else {
    return arb_pi();
  }
#endif
}

template <typename T> T real_euler() {
  if constexpr (std::is_same_v<T, double>) {
    return 0.57721566490153286061;
  } else if constexpr (std::is_same_v<T, DoubleDouble>) {
    return DoubleDoubleFromMpreal(mpfr::const_euler(128));
  } else if constexpr (std::is_same_v<T, mpfr::mpreal>) {
    return mpfr::const_euler();
  }
#ifdef ENABLE_FLINT_SUPPORT
  else {
    return arb_euler();
  }
#endif
}

// zeta(k), k >= 2
template <typename T> T real_zeta(int k) {
  if constexpr (std::is_same_v<T, double>) {
    return mpfr::zeta(mpfr::mpreal(k, 64)).toDouble();
  } else if constexpr (std::is_same_v<T, DoubleDouble>) {
    return DoubleDoubleFromMpreal(mpfr::zeta(mpfr::mpreal(k, 128)));
  } else if constexpr (std::is_same_v<T, mpfr::mpreal>) {
    return mpfr::zeta(mpfr::mpreal(k));
  }
#ifdef ENABLE_FLINT_SUPPORT
  else {
    return arb_zeta(static_cast<ulong>(k));
  }
#endif
}

template <typename T> T real_infinity() {
  if constexpr (std::is_same_v<T, double>) {
    return std::numeric_limits<double>::infinity();
  } else if constexpr (std::is_same_v<T, DoubleDouble>) {
    return DoubleDouble(std::numeric_limits<double>::infinity());
  } else if constexpr (std::is_same_v<T, mpfr::mpreal>) {
    return mpfr::const_infinity();
  }
#ifdef ENABLE_FLINT_SUPPORT
  else {
    ArbReal r;
    arb_pos_inf(r.v);
    return r;
  }
#endif
}

// x^e for an integer e, by binary exponentiation
template <typename T> T real_pow_int(T const &x, long e) {
  if (e < 0) {
    return real_from_int<T>(1) / real_pow_int(x, -e);
  }
  T res = real_from_int<T>(1);
  T b = x;
  while (e > 0) {
    if (e & 1) {
      res = res * b;
    }
    b = b * b;
    e >>= 1;
  }
  return res;
}

/*
  x in the type Tto, at the working precision of Tto: through the double,
  the DoubleDouble or the mpreal of x, the exact conversions being kept
  (a ball only goes to a ball, or to the floating types by its midpoint).
 */
template <typename Tto, typename Tfrom> Tto real_convert(Tfrom const &x) {
  if constexpr (std::is_same_v<Tto, Tfrom> &&
                !std::is_same_v<Tto, mpfr::mpreal>) {
    return x;
  } else if constexpr (std::is_same_v<Tto, double>) {
    return real_to_double(x);
  } else if constexpr (std::is_same_v<Tfrom, double>) {
    return real_from_double<Tto>(x);
  } else if constexpr (std::is_same_v<Tfrom, mpfr::mpreal>) {
    if constexpr (std::is_same_v<Tto, mpfr::mpreal>) {
      mpfr::mpreal y = x;
      y.setPrecision(static_cast<int>(mpfr::mpreal::get_default_prec()),
                     MPFR_RNDN);
      return y;
    } else if constexpr (std::is_same_v<Tto, DoubleDouble>) {
      return DoubleDoubleFromMpreal(x);
    }
#ifdef ENABLE_FLINT_SUPPORT
    else {
      // exact: the mpfr number as the midpoint
      ArbReal r;
      arf_set_mpfr(arb_midref(r.v), x.mpfr_srcptr());
      mag_zero(arb_radref(r.v));
      return r;
    }
#endif
  } else if constexpr (std::is_same_v<Tfrom, DoubleDouble>) {
    if constexpr (std::is_same_v<Tto, mpfr::mpreal>) {
      return DoubleDoubleToMpreal(x, mpfr::mpreal::get_default_prec());
    } else {
      return real_from_double<Tto>(x.hi) + real_from_double<Tto>(x.lo);
    }
  }
#ifdef ENABLE_FLINT_SUPPORT
  else {
    // from a ball, its midpoint
    mpfr::mpreal y(0, mpfr::mpreal::get_default_prec());
    arf_get_mpfr(y.mpfr_ptr(), arb_midref(x.v), MPFR_RNDN);
    if constexpr (std::is_same_v<Tto, DoubleDouble>) {
      return DoubleDoubleFromMpreal(y);
    } else {
      return y;
    }
  }
#endif
}

// clang-format off
#endif  // SRC_NUMBER_REALFUNCTIONS_H_
// clang-format on
