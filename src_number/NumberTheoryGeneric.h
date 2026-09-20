// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYGENERIC_H_
#define SRC_NUMBER_NUMBERTHEORYGENERIC_H_
// clang-format off
#include "Basic_functions.h"
#include "TemplateTraits.h"
#include <vector>
#include <optional>
#include <limits>
// clang-format on

template <typename T> T GenericGcd(T const &m, T const &n) {
  T h, q;
  if (n == 0 && m == 0) {
    return T(0);
  }
  T f = T_abs(m);
  T g = T_abs(n);
  while (g != 0) {
    q = QuoInt(f, g);
    h = g;
    g = f - q * g;
    f = h;
  }
  return f;
}

/*
  The inverse of a modulo P, or zero when a is not invertible.

  The classical extended Euclid carries a Bezout coefficient whose sign
  alternates, which is fine over the integers but wraps around on an
  unsigned type, where a negative value is not representable and the
  comparison against zero can never hold. So the recursion here is run on
  the magnitudes only,
     t_{k+1} = t_{k-1} + q_k t_k,
  which is non decreasing and stays bounded by P, so it neither goes
  negative nor leaves the range of a type that holds P. The sign that was
  dropped is recovered at the end from the number of steps, the identity
  being a t_k = (-1)^k gcd modulo P.

  The operands of the division are non negative throughout, so the
  truncation is the floor and the quotient is the Euclidean one for the
  signed and the unsigned types alike.
 */
template <typename T> T mod_inv(T const &a, T const &P) {
  T r0 = P;
  T r1 = ResInt(a, P);
  T t0(0);
  T t1(1);
  size_t n_step = 0;
  while (r1 != 0) {
    T q = r0 / r1;
    T r2 = r0 - q * r1;
    T t2 = t0 + q * t1;
    r0 = r1;
    r1 = r2;
    t0 = t1;
    t1 = t2;
    n_step++;
  }
  if (r0 != 1) {
    return T(0);
  }
  if (n_step % 2 == 1) {
    return t0;
  }
  return P - t0;
}

template <typename T>
requires (!is_mpz_class<T>::value)
inline PairGCD_dot<T> ComputePairGcdDot(T const &m, T const &n) {
  static_assert(is_euclidean_domain<T>::value,
                "Requires T to be an Euclidean domain in ComputePairGcd");
  T f, g, h, fm, gm, hm, q;
  if (n == 0 && m == 0) {
    f = 0;
    T a(0);
    T b(0);
    return {a, b, f};
  }
  if (m >= 0) {
    f = m;
    fm = 1;
  } else {
    f = -m;
    fm = -1;
  }
  if (n >= 0) {
    g = n;
    gm = 0;
  } else {
    g = -n;
    gm = 0;
  }
  while (g != 0) {
    q = QuoInt(f, g);
    h = g;
    hm = gm;
    g = f - q * g;
    gm = fm - q * gm;
    f = h;
    fm = hm;
  }
  T eCoeff1, eCoeff2;
  if (n == 0) {
    eCoeff1 = fm;
    eCoeff2 = 0;
  } else {
    eCoeff1 = fm;
    eCoeff2 = (f - fm * m) / n;
  }
#ifdef SANITY_CHECK_NB_THEORY_GENERIC
  T diff1 = f - eCoeff1 * m - eCoeff2 * n;
  if (diff1 != 0) {
    std::cerr << "A: diff1=" << diff1 << "\n";
    throw TerminalException{1};
  }
  // Check that coefficients are "small enough"
  // For Extended Euclidean Algorithm, coefficients should satisfy |eCoeff1| <=
  // |n/gcd| and |eCoeff2| <= |m/gcd| when both inputs are non-zero
  if (m != 0 && n != 0 && f != 0) {
    T abs_m = T_abs(m);
    T abs_n = T_abs(n);
    T abs_eCoeff1 = T_abs(eCoeff1);
    T abs_eCoeff2 = T_abs(eCoeff2);
    T abs_gcd = T_abs(f);

    T bound1 = abs_n / abs_gcd;
    T bound2 = abs_m / abs_gcd;

    if (abs_eCoeff1 > bound1) {
      std::cerr << "ERROR: |eCoeff1| = " << abs_eCoeff1 << " > " << bound1
                << " = |n|/|gcd|\n";
      std::cerr << "m=" << m << ", n=" << n << ", gcd=" << f << "\n";
      throw TerminalException{1};
    }

    if (abs_eCoeff2 > bound2) {
      std::cerr << "ERROR: |eCoeff2| = " << abs_eCoeff2 << " > " << bound2
                << " = |m|/|gcd|\n";
      std::cerr << "m=" << m << ", n=" << n << ", gcd=" << f << "\n";
      throw TerminalException{1};
    }
  }
#endif
  return {eCoeff1, eCoeff2, f};
}

template <typename T>
requires (!is_mpz_class<T>::value)
inline T KernelGcdPair(T const &a, T const &b) {
  return GenericGcd(a, b);
}

template <typename T>
requires is_totally_ordered<T>::value
inline T GcdPair(T const &a, T const &b) {
  T eGCD = KernelGcdPair(a, b);
  if (eGCD > 0)
    return eGCD;
  return -eGCD;
}

template <typename T>
requires (!is_totally_ordered<T>::value)
inline T GcdPair(T const &a, T const &b) {
  return KernelGcdPair(a, b);
}

template <typename T>
requires (!is_mpz_class<T>::value)
inline T KernelLCMpair(T const &a, T const &b) {
  if (a == 0)
    return b;
  if (b == 0)
    return a;
  return a * b / KernelGcdPair(a, b);
}

template <typename T>
requires (!is_totally_ordered<T>::value)
inline T LCMpair(T const &a, T const &b) {
  return KernelLCMpair(a, b);
}

template <typename T>
requires is_totally_ordered<T>::value
inline T LCMpair(T const &a, T const &b) {
  T eLCM = KernelLCMpair(a, b);
  if (eLCM > 0)
    return eLCM;
  return -eLCM;
}

template <typename T> T LCMlist(std::vector<T> const &V) {
  size_t len = V.size();
  T eLCM = V[0];
  for (size_t u = 1; u < len; u++) {
    eLCM = LCMpair(eLCM, V[u]);
  }
  return eLCM;
}

template <typename T> std::optional<T> UniversalSquareRoot(T const &val) {
  if (val < 0)
    return {};
  T ret;
  if (!universal_square_root(ret, val))
    return {};
  return ret;
}

/*
  Given a vector of a and a vector of m find a x such that
  x = a[i] mod m[i] for all i
  the m[i] need to be coprime.
  ---
  We apply
  https://en.wikipedia.org/wiki/Chinese_remainder_theorem
 */
template <typename T>
T chinese_remainder_theorem(std::vector<T> const &a, std::vector<T> const &m) {
#ifdef DEBUG_NUMBER_THEORY_GENERIC
  if (a.size() != m.size()) {
    std::cerr << "a and m should be of equal lengths\n";
    throw TerminalException{1};
  }
  if (a.size() == 0) {
    std::cerr << "a should be of positive length\n";
    throw TerminalException{1};
  }
#endif
  size_t siz = m.size();
  T x = a[0];
  T m_prod = m[0];
  for (size_t i = 1; i < siz; i++) {
    PairGCD_dot<T> t = ComputePairGcdDot(m_prod, m[i]);
#ifdef DEBUG_NUMBER_THEORY_GENERIC
    if (t.gcd != 1) {
      std::cerr << "The GCD should be equal to 1\n";
      throw TerminalException{1};
    }
#endif
    x = x * t.b * m[i] + a[i] * t.a * m_prod;
#ifdef DEBUG_NUMBER_THEORY_GENERIC
    T sum = t.b * m[i] + t.a * m_prod;
    if (sum != 1) {
      std::cerr << "The t is not correct\n";
      throw TerminalException{1};
    }
#endif
    m_prod *= m[i];
  }
#ifdef DEBUG_NUMBER_THEORY_GENERIC
  for (size_t i = 0; i < siz; i++) {
    T diff = x - a[i];
    T res = ResInt(diff, m[i]);
    if (res != 0) {
      std::cerr << "NTG: a=";
      for (auto &val : a) {
        std::cerr << " " << val;
      }
      std::cerr << "\n";
      std::cerr << "NTG: m=";
      for (auto &val : m) {
        std::cerr << " " << val;
      }
      std::cerr << "\n";
      std::cerr << "NTG: x=" << x << "\n";
      std::cerr << "NTG: We do not have a solution of the Chinese Remainder "
                   "Theorem\n";
      throw TerminalException{1};
    }
  }
#endif
  return x;
}

template <typename T>
requires std::is_integral<T>::value
inline void set_to_infinity(T &x) {
  x = std::numeric_limits<T>::max();
}

template <typename T> T practical_infinity() {
  T ret;
  set_to_infinity(ret);
  return ret;
}

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYGENERIC_H_
// clang-format on
