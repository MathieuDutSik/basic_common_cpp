// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_ARITHMETICFUNCTIONS_H_
#define SRC_NUMBER_ARITHMETICFUNCTIONS_H_

// clang-format off
#include "factorizations.h"
#include <map>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

/*
  The arithmetic functions of the integers read off their factorization
  (FactorsIntMap, factorizations.h), for any integer type T:
  - EulerPhi(n), Moebius(n) for n >= 1;
  - IsSquareFree(n), SquareFreePart(n) (with the sign of n) for n != 0;
  - PrimeDivisors(n), OddPrimeDivisors(n): the primes dividing n != 0,
    without multiplicity and sorted;
  - PrimePowerDecomposition(q): the prime p and the exponent k >= 1 with
    q = p^k, nothing when q is not such a power;
  - BinomialCoefficient<T>(n, k), computed in T from integers n, k of
    another type, and the Gaussian binomial [a b]_q.
  The factorization is the expensive part, so a caller needing several of
  them for the same n can factor once and read the map.
 */

template <typename T> T IntegerAbsoluteValue(T const &n) {
  return (n < 0) ? T(-n) : n;
}

template <typename T> T EulerPhi(T const &n) {
  T res(1);
  if (n <= 1) {
    return res;
  }
  for (auto const &[p, mult] : FactorsIntMap(n)) {
    res *= p - 1;
    for (size_t i = 1; i < mult; i++) {
      res *= p;
    }
  }
  return res;
}

template <typename T> int Moebius(T const &n) {
  if (n <= 1) {
    return 1;
  }
  int mu = 1;
  for (auto const &[p, mult] : FactorsIntMap(n)) {
    if (mult > 1) {
      return 0;
    }
    mu = -mu;
  }
  return mu;
}

template <typename T> bool IsSquareFree(T const &n) {
  T n_abs = IntegerAbsoluteValue(n);
  if (n_abs <= 1) {
    return true;
  }
  for (auto const &[p, mult] : FactorsIntMap(n_abs)) {
    if (mult > 1) {
      return false;
    }
  }
  return true;
}

// The square-free part of n != 0, with the sign of n.
template <typename T> T SquareFreePart(T const &n) {
  T res(1);
  T n_abs = IntegerAbsoluteValue(n);
  if (n_abs > 1) {
    for (auto const &[p, mult] : FactorsIntMap(n_abs)) {
      if (mult % 2 == 1) {
        res *= p;
      }
    }
  }
  return (n < 0) ? T(-res) : res;
}

template <typename T> std::vector<T> PrimeDivisors(T const &n) {
  std::vector<T> res;
  T n_abs = IntegerAbsoluteValue(n);
  if (n_abs > 1) {
    for (auto const &[p, mult] : FactorsIntMap(n_abs)) {
      res.push_back(p);
    }
  }
  return res;
}

template <typename T> std::vector<T> OddPrimeDivisors(T const &n) {
  std::vector<T> res;
  for (auto const &p : PrimeDivisors(n)) {
    if (p != 2) {
      res.push_back(p);
    }
  }
  return res;
}

template <typename T>
std::optional<std::pair<T, size_t>> PrimePowerDecomposition(T const &q) {
  if (q <= 1) {
    return {};
  }
  std::map<T, size_t> map = FactorsIntMap(q);
  if (map.size() != 1) {
    return {};
  }
  return *map.begin();
}

// C(n, k) = (n / k) C(n - 1, k - 1), 0 for k < 0.
template <typename T, typename Tint>
T BinomialCoefficient(Tint const &n, Tint const &k) {
  T res(1);
  if (k < 0) {
    return T(0);
  }
  for (Tint i = 1; i <= k; i++) {
    T num = n - (k - i);
    T den = i;
    res *= num;
    res /= den;
  }
  return res;
}

/*
  The Gaussian binomial [a b]_q = prod_{i=1}^{b} (q^(a-b+i) - 1) / (q^i - 1),
  0 unless 0 <= b <= a, by the symmetry [a b] = [a a-b] with b <= a / 2:
  the powers of q are carried along and the quotient is taken once.
 */
template <typename T> T GaussianBinomial(int a, int b, T const &q) {
  if (b < 0 || b > a) {
    return T(0);
  }
  if (2 * b > a) {
    b = a - b;
  }
  T q_num(1);
  for (int i = 0; i < a - b; i++) {
    q_num *= q;
  }
  T q_den(1), num(1), den(1);
  for (int i = 1; i <= b; i++) {
    q_num *= q;
    q_den *= q;
    num *= q_num - 1;
    den *= q_den - 1;
  }
  return num / den;
}

// clang-format off
#endif  // SRC_NUMBER_ARITHMETICFUNCTIONS_H_
// clang-format on
