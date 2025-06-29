// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYGENERIC_H_
#define SRC_NUMBER_NUMBERTHEORYGENERIC_H_
// clang-format off
#include "Basic_functions.h"
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

template <typename T> T mod_inv(T const &a, T const &P) {
  T t(0);
  T newt(1);
  T r = P;
  T newr = a;
  T q, tmp;
  while (newr != 0) {
    q = QuoInt(r, newr);
    tmp = t;
    t = newt;
    newt = tmp - q * newt;
    tmp = r;
    r = newr;
    newr = tmp - q * newr;
  }
  if (r > 1)
    return T(0);
  if (t < 0)
    t = t + P;
  return t;
}

template <typename T>
inline typename std::enable_if<!is_mpz_class<T>::value, PairGCD_dot<T>>::type
ComputePairGcdDot(T const &m, T const &n) {
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
#ifdef DEBUG_MATRIX_INT
  T diff1 = f - eCoeff1 * m - eCoeff1 * n;
  if (diff1 != 0) {
    std::cerr << "A: diff1=" << diff1 << "\n";
    throw TerminalException{1};
  }
#endif
  return {eCoeff1, eCoeff2, f};
}

template <typename T>
inline typename std::enable_if<!is_mpz_class<T>::value, T>::type
KernelGcdPair(T const &a, T const &b) {
  return GenericGcd(a, b);
}

template <typename T>
inline typename std::enable_if<is_totally_ordered<T>::value, T>::type
GcdPair(T const &a, T const &b) {
  T eGCD = KernelGcdPair(a, b);
  if (eGCD > 0)
    return eGCD;
  return -eGCD;
}

template <typename T>
inline typename std::enable_if<!is_totally_ordered<T>::value, T>::type
GcdPair(T const &a, T const &b) {
  return KernelGcdPair(a, b);
}

template <typename T>
inline typename std::enable_if<!is_mpz_class<T>::value, T>::type
KernelLCMpair(T const &a, T const &b) {
  if (a == 0)
    return b;
  if (b == 0)
    return a;
  return a * b / KernelGcdPair(a, b);
}

template <typename T>
inline typename std::enable_if<!is_totally_ordered<T>::value, T>::type
LCMpair(T const &a, T const &b) {
  return KernelLCMpair(a, b);
}

template <typename T>
inline typename std::enable_if<is_totally_ordered<T>::value, T>::type
LCMpair(T const &a, T const &b) {
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
template<typename T>
T chinese_remainder_theorem(std::vector<T> const& a, std::vector<T> const& m) {
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
  for (size_t i=1; i<siz; i++) {
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
  for (size_t i=0; i<siz; i++) {
    T diff = x - a[i];
    T res = ResInt(diff, m[i]);
    if (res != 0) {
      std::cerr << "NTG: a=";
      for (auto & val : a) {
        std::cerr << " " << val;
      }
      std::cerr << "\n";
      std::cerr << "NTG: m=";
      for (auto & val : m) {
        std::cerr << " " << val;
      }
      std::cerr << "\n";
      std::cerr << "NTG: x=" << x << "\n";
      std::cerr << "NTG: We do not have a solution of the Chinese Remainder Theorem\n";
      throw TerminalException{1};
    }
  }
#endif
  return x;
}



template <typename T>
inline typename std::enable_if<std::is_integral<T>::value, void>::type
set_to_infinity(T &x) {
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
