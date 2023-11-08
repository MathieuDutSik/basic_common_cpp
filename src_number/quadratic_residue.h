// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_FACTORIZATIONS_H_
#define SRC_NUMBER_FACTORIZATIONS_H_

// clang-format off
#include "TemplateTraits.h"
#include <map>
#include <utility>
#include <vector>
// clang-format on

template <typename T>
std::optional<T> find_quadratic_residue(T const& a, T const& m_in) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  T m = T_abs(m_in);
  T two(2);
  T a_mod = QuoInt(a, m);
  T res = ResInt(m, two);
  T upper(0);
  if (res == 0) {
    upper = QuoInt(m, two);
  } else {
    T mP1 = m+1;
    upper = QuoInt(mP1, two);
  }
  T x(0);
  T TwoXpOne(1);
  T xSqr(0);
  while (x != upper) {
    if (xSqr == a_mod) {
      return x;
    }
    xSqr += TwoXpOne;
    xSqr = QuoInt(xSqr, m);
    TwoXpOne += 2;
    x += 2;
  }
  return {};
}

template <typename T>
bool is_quadratic_residue(T const& a, T const& m) {
  std::optional<T> opt = find_quadratic_residue(a, m);
  return opt.has_value();
}

// clang-format off
#endif  // SRC_NUMBER_FACTORIZATIONS_H_
// clang-format on
