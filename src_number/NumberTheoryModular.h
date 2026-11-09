// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYMODULAR_H_
#define SRC_NUMBER_NUMBERTHEORYMODULAR_H_
// clang-format off
#include "NumberTheoryGeneric.h"
#include "QuoIntFcts.h"
#include <optional>
#include <vector>
// clang-format on

/*
  Modular and p-adic operations on the integer types (int64_t, mpz_class,
  fmpz_class, the boost integers, and the rational types used for integer
  values), written with ResInt and QuoInt only, so that they run on every
  type that has them, with no use of a library-specific call.

  - power_modulo(a, e, m): a^e modulo m, e >= 0, m > 0, in [0, m).
  - jacobi_symbol(a, n): the Jacobi symbol (a / n) for n odd and positive,
    0 when gcd(a, n) > 1.
  - kronecker_symbol(a, n): the Kronecker symbol (a / n) for n > 0.
  - prime_valuation(x, p): v_p(x) for x != 0 and p > 1.
  - rational_prime_valuation(x, p): v_p of a nonzero rational, its
    numerator and denominator being read with GetNumerator_z and
    GetDenominator_z.
  - exact_integer_root(x, r): the y with y^r = x when there is one (x >= 0,
    or r odd), by Newton's iteration on the integers.
 */

template <typename T> T power_modulo(T const &a, T const &e, T const &m) {
  T result = ResInt(T(1), m);
  T base = ResInt(a, m);
  T exp = e;
  T two(2);
  while (exp > 0) {
    if (ResInt(exp, two) == 1) {
      T prod = result * base;
      result = ResInt(prod, m);
    }
    T sq = base * base;
    base = ResInt(sq, m);
    exp = QuoInt(exp, two);
  }
  return result;
}

template <typename T> int jacobi_symbol(T const &a_in, T const &n_in) {
  T two(2), four(4), eight(8);
  T a = ResInt(a_in, n_in);
  T n = n_in;
  int result = 1;
  while (a != 0) {
    while (ResInt(a, two) == 0) {
      a = QuoInt(a, two);
      T r = ResInt(n, eight);
      if (r == 3 || r == 5) {
        result = -result;
      }
    }
    T t = a;
    a = n;
    n = t;
    if (ResInt(a, four) == 3 && ResInt(n, four) == 3) {
      result = -result;
    }
    a = ResInt(a, n);
  }
  return (n == 1) ? result : 0;
}

template <typename T> int kronecker_symbol(T const &a, T const &n_in) {
  T two(2), eight(8);
  T n = n_in;
  int result = 1;
  while (ResInt(n, two) == 0) {
    n = QuoInt(n, two);
    T r = ResInt(a, eight);
    if (ResInt(r, two) == 0) {
      return 0;
    }
    if (r == 3 || r == 5) {
      result = -result;
    }
  }
  if (n == 1) {
    return result;
  }
  return result * jacobi_symbol(a, n);
}

template <typename T> size_t prime_valuation(T const &x, T const &p) {
  size_t v = 0;
  T y = x;
  while (y != 0 && ResInt(y, p) == 0) {
    y = QuoInt(y, p);
    v++;
  }
  return v;
}

template <typename Tq, typename Tz>
int rational_prime_valuation(Tq const &x, Tz const &p) {
  Tz num = GetNumerator_z(x);
  Tz den = GetDenominator_z(x);
  return static_cast<int>(prime_valuation(num, p)) -
         static_cast<int>(prime_valuation(den, p));
}

template <typename T>
std::optional<T> exact_integer_root(T const &x, size_t r) {
  if (r == 0) {
    return {};
  }
  if (x < 0) {
    if (r % 2 == 0) {
      return {};
    }
    T mx = -x;
    std::optional<T> y = exact_integer_root(mx, r);
    if (!y) {
      return {};
    }
    T my = -*y;
    return my;
  }
  if (x == 0 || x == 1 || r == 1) {
    return x;
  }
  auto power = [&](T const &y) -> T {
    T res(1);
    for (size_t i = 0; i < r; i++) {
      res *= y;
    }
    return res;
  };
  // Newton from above: y <- ((r - 1) y + x / y^(r-1)) / r, decreasing to
  // the floor of the root
  T y(1);
  while (power(y) < x) {
    y *= 2;
  }
  T r_T(0), rm1_T(0);
  for (size_t i = 0; i < r; i++) {
    r_T += 1;
  }
  rm1_T = r_T - 1;
  while (true) {
    T ym1(1);
    for (size_t i = 0; i + 1 < r; i++) {
      ym1 *= y;
    }
    T num = rm1_T * y;
    num += QuoInt(x, ym1);
    T z = QuoInt(num, r_T);
    if (z >= y) {
      break;
    }
    y = z;
  }
  if (power(y) == x) {
    return y;
  }
  return {};
}

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYMODULAR_H_
// clang-format on
