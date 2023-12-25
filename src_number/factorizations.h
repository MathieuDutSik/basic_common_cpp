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
std::pair<bool, T> rho_pollard_factorize(T const &number) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  T count;
  T x_fixed = 2, x = 2, size = 2, factor, diff;
  do {
    count = size;
    do {
      ResInt_Kernel(x * x + 1, number, x);
      diff = x - x_fixed;
      if (diff < 0)
        diff = -diff;
      factor = GcdPair(diff, number);
    } while (--count > 0 && factor == 1);
    size *= 2;
    x_fixed = x;
  } while (factor == 1);
  if (factor == number) {
    return {false, -1};
  } else {
    return {true, factor};
  }
}

template <typename T> std::vector<T> successive_division_factorize(T const &N) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  T pos = 2;
  while (true) {
    T res = ResInt(N, pos);
    if (res == 0) {
      T quot = QuoInt(N, pos);
      if (quot > 1) {
        std::vector<T> eVect = successive_division_factorize(quot);
        eVect.push_back(pos);
        return eVect;
      }
      return {pos};
    }
    pos++;
    if (pos * pos > N)
      break;
  }
  return {N};
}

template <typename T> bool successive_division_isprime(T const &N) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  T pos = 2;
  while (true) {
    T res = ResInt(N, pos);
    if (res == 0)
      return false;
    pos++;
    if (pos * pos > N)
      break;
  }
  return true;
}

template <typename T> bool IsPrime(const T &N) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  std::pair<bool, T> epair = rho_pollard_factorize(N);
  if (epair.first) {
    return false;
  } else {
    return successive_division_isprime(N);
  }
}

template <typename T> std::vector<T> Kernel_FactorsInt(T const &N) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  if (N == 1)
    return {};
  std::pair<bool, T> epair = rho_pollard_factorize(N);
  if (epair.first) {
    T fact1 = epair.second;
    T fact2 = QuoInt(N, fact1);
    std::vector<T> ListPrime = Kernel_FactorsInt(fact1);
    std::vector<T> V2 = Kernel_FactorsInt(fact2);
    ListPrime.insert(ListPrime.end(), V2.begin(), V2.end());
    return ListPrime;
  } else {
    return successive_division_factorize(N);
  }
}

template <typename T>
inline typename std::enable_if<is_implementation_of_Z<T>::value,
                               std::vector<T>>::type
FactorsInt(T const &N) {
  return Kernel_FactorsInt(N);
}

template <typename T>
inline typename std::enable_if<!is_implementation_of_Z<T>::value,
                               std::vector<T>>::type
FactorsInt(T const &N) {
  using Tint = typename underlying_ring<T>::ring_type;
  Tint N_int = UniversalScalarConversion<Tint, T>(N);
  std::vector<Tint> LFact_int = Kernel_FactorsInt(N_int);
  std::vector<T> LFact;
  for (auto &val_i : LFact_int) {
    T val = UniversalScalarConversion<T, Tint>(val_i);
    LFact.push_back(val);
  }
  return LFact;
}

template <typename T> std::map<T, size_t> FactorsIntMap(T const &N) {
  std::vector<T> vect = FactorsInt(N);
  std::map<T, size_t> map;
  for (auto &eV : vect) {
    map[eV] += 1;
  }
  return map;
}

template <typename T>
std::vector<T> GetAllFactors(std::map<T, int> const &eMap) {
  std::vector<T> LVal = {1};
  for (auto &kv : eMap) {
    std::vector<T> NewVal;
    T ePow = 1;
    T mult = kv.first;
    for (int i = 0; i <= kv.second; i++) {
      for (auto &eVal : LVal)
        NewVal.push_back(ePow * eVal);
      ePow *= mult;
    }
    LVal = NewVal;
  }
  return LVal;
}

template <typename T> std::vector<T> GetAllFactors(T const &N) {
  std::vector<T> LFact = FactorsInt(N);
  std::map<T, int> eMap;
  for (auto &eVal : LFact)
    eMap[eVal]++;
  return GetAllFactors(eMap);
}

// clang-format off
#endif  // SRC_NUMBER_FACTORIZATIONS_H_
// clang-format on
