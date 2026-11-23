// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_DOUBLEDOUBLE_H_
#define SRC_NUMBER_DOUBLEDOUBLE_H_

// clang-format off
#include <cmath>
#include <cstdint>
// clang-format on

/*
  Double-double numbers hi + lo, |lo| <= ulp(hi) / 2, about 106 bits
  (Dekker, Knuth; Hida, Li and Bailey): the field operations, exp and
  log. Faster than a multiprecision type for the computations needing a
  little more than double precision.
 */
struct DoubleDouble {
  double hi;
  double lo;

  DoubleDouble() : hi(0.0), lo(0.0) {}
  DoubleDouble(double h) : hi(h), lo(0.0) {}
  DoubleDouble(double h, double l) : hi(h), lo(l) {}
  explicit DoubleDouble(int64_t x) {
    // exact: x = hi + (x - hi), the remainder below 2^11
    hi = static_cast<double>(x);
    lo = static_cast<double>(x - static_cast<int64_t>(hi));
  }
  explicit DoubleDouble(int x) : hi(x), lo(0.0) {}
};

inline DoubleDouble DoubleDoubleTwoSum(double a, double b) {
  double s = a + b;
  double bb = s - a;
  double e = (a - (s - bb)) + (b - bb);
  return {s, e};
}

inline DoubleDouble DoubleDoubleQuickTwoSum(double a, double b) {
  double s = a + b;
  return {s, b - (s - a)};
}

inline DoubleDouble operator+(DoubleDouble const &a, DoubleDouble const &b) {
  DoubleDouble s = DoubleDoubleTwoSum(a.hi, b.hi);
  DoubleDouble t = DoubleDoubleTwoSum(a.lo, b.lo);
  s.lo += t.hi;
  s = DoubleDoubleQuickTwoSum(s.hi, s.lo);
  s.lo += t.lo;
  return DoubleDoubleQuickTwoSum(s.hi, s.lo);
}

inline DoubleDouble operator-(DoubleDouble const &a) { return {-a.hi, -a.lo}; }

inline DoubleDouble operator-(DoubleDouble const &a, DoubleDouble const &b) {
  return a + (-b);
}

inline DoubleDouble operator*(DoubleDouble const &a, DoubleDouble const &b) {
  double p = a.hi * b.hi;
  double e = std::fma(a.hi, b.hi, -p);
  e += a.hi * b.lo + a.lo * b.hi;
  return DoubleDoubleQuickTwoSum(p, e);
}

// a / b: q1 = a / b in double, then two corrections by the remainder
inline DoubleDouble operator/(DoubleDouble const &a, DoubleDouble const &b) {
  double q1 = a.hi / b.hi;
  DoubleDouble r = a - b * DoubleDouble(q1);
  double q2 = r.hi / b.hi;
  r = r - b * DoubleDouble(q2);
  double q3 = r.hi / b.hi;
  DoubleDouble q = DoubleDoubleQuickTwoSum(q1, q2);
  return q + DoubleDouble(q3);
}

inline DoubleDouble &operator+=(DoubleDouble &a, DoubleDouble const &b) {
  a = a + b;
  return a;
}
inline DoubleDouble &operator-=(DoubleDouble &a, DoubleDouble const &b) {
  a = a - b;
  return a;
}
inline DoubleDouble &operator*=(DoubleDouble &a, DoubleDouble const &b) {
  a = a * b;
  return a;
}
inline DoubleDouble &operator/=(DoubleDouble &a, DoubleDouble const &b) {
  a = a / b;
  return a;
}

inline bool operator<(DoubleDouble const &a, DoubleDouble const &b) {
  return a.hi < b.hi || (a.hi == b.hi && a.lo < b.lo);
}
inline bool operator>(DoubleDouble const &a, DoubleDouble const &b) {
  return b < a;
}
inline bool operator<=(DoubleDouble const &a, DoubleDouble const &b) {
  return !(b < a);
}
inline bool operator>=(DoubleDouble const &a, DoubleDouble const &b) {
  return !(a < b);
}
inline bool operator==(DoubleDouble const &a, DoubleDouble const &b) {
  return a.hi == b.hi && a.lo == b.lo;
}

// exp(x) for |x| < 700: x = k log 2 + r, exp(r / 2^10) by Taylor, then 10
// squarings kept as expm1 (2 e + e^2), and the scaling by 2^k
inline DoubleDouble DoubleDoubleExp(DoubleDouble const &x) {
  const DoubleDouble ln2{0.6931471805599453, 2.3190468138462996e-17};
  double k = std::nearbyint(x.hi / ln2.hi);
  DoubleDouble r = x - ln2 * DoubleDouble(k);
  r.hi = std::ldexp(r.hi, -10);
  r.lo = std::ldexp(r.lo, -10);
  // expm1(r) = r + r^2/2 + ... + r^9/9!
  DoubleDouble e = r;
  DoubleDouble term = r;
  for (int i = 2; i <= 9; i++) {
    term = term * r;
    term.hi /= i;
    term.lo /= i;
    term = DoubleDoubleQuickTwoSum(term.hi, term.lo);
    e = e + term;
  }
  for (int i = 0; i < 10; i++) {
    // expm1(2y) = 2 expm1(y) + expm1(y)^2
    e = DoubleDouble(2 * e.hi, 2 * e.lo) + e * e;
  }
  DoubleDouble res = DoubleDouble(1.0) + e;
  int ki = static_cast<int>(k);
  return {std::ldexp(res.hi, ki), std::ldexp(res.lo, ki)};
}

// log x for x > 0: y = log x in double, then one Newton step y + x e^{-y}
// - 1, which doubles the 53 bits
inline DoubleDouble DoubleDoubleLog(DoubleDouble const &x) {
  double y = std::log(x.hi);
  DoubleDouble t = x * DoubleDoubleExp(DoubleDouble(-y));
  return DoubleDouble(y) + (t - DoubleDouble(1.0));
}

// clang-format off
#endif  // SRC_NUMBER_DOUBLEDOUBLE_H_
// clang-format on
