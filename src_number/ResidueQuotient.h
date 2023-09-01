// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_RESIDUEQUOTIENT_H_
#define SRC_NUMBER_RESIDUEQUOTIENT_H_

#include "BasicNumberTypes.h"
#include <cstdint>
#include <cstdlib>

// The remainder and quotient of integers

template <typename T> T ResInt_C_integer(T const &a, T const &b) {
  T res2 = a % b;
  if (a < 0 && res2 != 0)
    res2 += std::abs(b);
  return res2;
}

template <typename T> T ResInt_C_unsigned_integer(T const &a, T const &b) {
  T res2 = a % b;
  return res2;
}

// Specific functions for number type.
// We basically want to do the PID (Principal Ideal Domain) case.
// What are needed for that are two functions:
// ---QuotInt for the quotient
// ---T_Norm for the norm of elements.
// Result of QuoInt(a,b) is an integer q such that q is the quotient.
// We then have a = q b + r
//
// For natural integer Z (i.e. int/long/mpz_class/mpq_class)
// We should have a = bq + r
// with 0 <= r < |b| and q integer.
template <typename T> T ResInt_Generic(T const &a, T const &b) {
  // We cannot use std::abs which is not defined for all data types.
  T b_abs;
  if (b > 0)
    b_abs = b;
  else
    b_abs = -b;
  T res = a % b_abs;
  while (true) {
    if (res >= 0 && res < b_abs)
      break;
    if (res < 0)
      res += b_abs;
    if (res >= b_abs)
      res -= b_abs;
  }
  return res;
}

inline void ResInt_Kernel(int const &a, int const &b, int &res) {
  res = ResInt_C_integer<int>(a, b);
}

inline void ResInt_Kernel(uint8_t const &a, uint8_t const &b, uint8_t &res) {
  res = ResInt_C_unsigned_integer<uint8_t>(a, b);
}

inline void ResInt_Kernel(uint16_t const &a, uint16_t const &b, uint16_t &res) {
  res = ResInt_C_unsigned_integer<uint16_t>(a, b);
}

inline void ResInt_Kernel(uint32_t const &a, uint32_t const &b, uint32_t &res) {
  res = ResInt_C_unsigned_integer<uint32_t>(a, b);
}

inline void ResInt_Kernel(long const &a, long const &b, long &res) {
  res = ResInt_C_integer<long>(a, b);
}

template <typename T> T QuoInt_C_integer(T const &a, T const &b) {
  T quo2 = a / b;
  if (a < 0 && b * quo2 != a) {
    if (b > 0)
      quo2--;
    else
      quo2++;
  }
  return quo2;
}

template <typename T> T QuoInt_Generic(T const &a, T const &b) {
  T res = ResInt_Generic(a, b);
  return (a - res) / b;
}

inline void QUO_INT(stc<int> const &a, stc<int> const &b, int & q) {
  q = QuoInt_C_integer<int>(a.val, b.val);
}

inline void QUO_INT(stc<long> const &a, stc<long> const &b, long & q) {
  q = QuoInt_C_integer<long>(a.val, b.val);
}

#include "QuoIntFcts.h"

inline int GetDenominator([[maybe_unused]] int const &x) { return 1; }

inline long GetDenominator([[maybe_unused]] long const &x) { return 1; }

inline int GetNumerator_z(int const &x) { return x; }

inline int GetDenominator_z([[maybe_unused]] int const &x) { return 1; }

inline long GetNumerator_z(long const &x) { return x; }

inline long GetDenominator_z([[maybe_unused]] long const &x) { return 1; }

template <typename T> T ResInt(T const &a, T const &b) {
  T res;
  ResInt_Kernel(a, b, res);
  return res;
}

// clang-format off
#endif  // SRC_NUMBER_RESIDUEQUOTIENT_H_
// clang-format on
