// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYGENERIC_H_
#define SRC_NUMBER_NUMBERTHEORYGENERIC_H_

#include "Basic_functions.h"
#include <vector>

template <typename T> T GenericGcd(T const &m, T const &n) {
  T h, q;
  if (n == 0 && m == 0) {
    return 0;
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
