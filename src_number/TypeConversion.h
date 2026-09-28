// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_TYPECONVERSION_H_
#define SRC_NUMBER_TYPECONVERSION_H_

// clang-format off
#include "BasicNumberTypes.h"
#include "ExceptionsFunc.h"
#include "TemplateTraits.h"
#include <cmath>
#include <cstdint>
#include <iostream>
#include <math.h>
#include <type_traits>
#include <utility>
#include <vector>
#include <optional>
// clang-format on

// All the definitions of special fields are in other include.
// Nothing of this should depend on GMP or MPREAL or FLINT or whatever.
//
//  All mpreal are in mpreal_related.h

//
// UniversalScalarConversion and TYPE_CONVERSION
//

// Conversion from double
inline void TYPE_CONVERSION(stc<double> const &a1, double &a2) { a2 = a1.val; }

inline void TYPE_CONVERSION(stc<double> const &a1, uint8_t &a2) {
  a2 = static_cast<uint8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<double> const &a1, int8_t &a2) {
  a2 = static_cast<int8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<double> const &a1, uint16_t &a2) {
  a2 = static_cast<uint16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<double> const &a1, int16_t &a2) {
  a2 = static_cast<int16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<double> const &a1, uint32_t &a2) {
  a2 = static_cast<uint32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<double> const &a1, int32_t &a2) {
  a2 = static_cast<int32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<double> const &a1, uint64_t &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<double> const &a1, int64_t &a2) {
  a2 = static_cast<int64_t>(a1.val);
}

// Conversion from int8_t

inline void TYPE_CONVERSION(stc<int8_t> const &a1, double &a2) {
  a2 = static_cast<double>(a1.val);
}

inline void TYPE_CONVERSION(stc<int8_t> const &a1, uint8_t &a2) {
  a2 = static_cast<uint8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int8_t> const &a1, int8_t &a2) {
  a2 = static_cast<int8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int8_t> const &a1, uint16_t &a2) {
  a2 = static_cast<uint16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int8_t> const &a1, int16_t &a2) {
  a2 = static_cast<int16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int8_t> const &a1, uint32_t &a2) {
  a2 = static_cast<uint32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int8_t> const &a1, int32_t &a2) {
  a2 = static_cast<int32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int8_t> const &a1, uint64_t &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<int8_t> const &a1, int64_t &a2) {
  a2 = static_cast<int64_t>(a1.val);
}

// Conversion from uint8_t

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, double &a2) {
  a2 = static_cast<double>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, uint8_t &a2) {
  a2 = static_cast<uint8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, int8_t &a2) {
  a2 = static_cast<int8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, uint16_t &a2) {
  a2 = static_cast<uint16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, int16_t &a2) {
  a2 = static_cast<int16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, uint32_t &a2) {
  a2 = static_cast<uint32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, int32_t &a2) {
  a2 = static_cast<int32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, uint64_t &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, int64_t &a2) {
  a2 = static_cast<int64_t>(a1.val);
}

// Conversion from int16_t

inline void TYPE_CONVERSION(stc<int16_t> const &a1, double &a2) {
  a2 = static_cast<double>(a1.val);
}

inline void TYPE_CONVERSION(stc<int16_t> const &a1, uint8_t &a2) {
  a2 = static_cast<uint8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int16_t> const &a1, int8_t &a2) {
  a2 = static_cast<int8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int16_t> const &a1, uint16_t &a2) {
  a2 = static_cast<uint16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int16_t> const &a1, int16_t &a2) {
  a2 = static_cast<int16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int16_t> const &a1, uint32_t &a2) {
  a2 = static_cast<uint32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int16_t> const &a1, int32_t &a2) {
  a2 = static_cast<int32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int16_t> const &a1, uint64_t &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<int16_t> const &a1, int64_t &a2) {
  a2 = static_cast<int64_t>(a1.val);
}

// Conversion from uint16_t

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, double &a2) {
  a2 = static_cast<double>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, uint8_t &a2) {
  a2 = static_cast<uint8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, int8_t &a2) {
  a2 = static_cast<int8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, uint16_t &a2) {
  a2 = static_cast<uint16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, int16_t &a2) {
  a2 = static_cast<int16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, uint32_t &a2) {
  a2 = static_cast<uint32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, int32_t &a2) {
  a2 = static_cast<int32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, uint64_t &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, int64_t &a2) {
  a2 = static_cast<int64_t>(a1.val);
}

// Conversion from int32_t

inline void TYPE_CONVERSION(stc<int32_t> const &a1, double &a2) {
  a2 = static_cast<double>(a1.val);
}

inline void TYPE_CONVERSION(stc<int32_t> const &a1, uint8_t &a2) {
  a2 = static_cast<uint8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int32_t> const &a1, int8_t &a2) {
  a2 = static_cast<int8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int32_t> const &a1, uint16_t &a2) {
  a2 = static_cast<uint16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int32_t> const &a1, int16_t &a2) {
  a2 = static_cast<int16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int32_t> const &a1, uint32_t &a2) {
  a2 = static_cast<uint32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int32_t> const &a1, int32_t &a2) {
  a2 = static_cast<int32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int32_t> const &a1, uint64_t &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<int32_t> const &a1, int64_t &a2) {
  a2 = static_cast<int64_t>(a1.val);
}

// Conversion from uint32_t

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, double &a2) {
  a2 = static_cast<double>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, uint8_t &a2) {
  a2 = static_cast<uint8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, int8_t &a2) {
  a2 = static_cast<int8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, uint16_t &a2) {
  a2 = static_cast<uint16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, int16_t &a2) {
  a2 = static_cast<int16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, uint32_t &a2) {
  a2 = static_cast<uint32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, int32_t &a2) {
  a2 = static_cast<int32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, uint64_t &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, int64_t &a2) {
  a2 = static_cast<int64_t>(a1.val);
}

// Conversion from int64_t

inline void TYPE_CONVERSION(stc<int64_t> const &a1, double &a2) {
  a2 = static_cast<double>(a1.val);
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, uint8_t &a2) {
  a2 = static_cast<uint8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, int8_t &a2) {
  a2 = static_cast<int8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, uint16_t &a2) {
  a2 = static_cast<uint16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, int16_t &a2) {
  a2 = static_cast<int16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, uint32_t &a2) {
  a2 = static_cast<uint32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, int32_t &a2) {
  a2 = static_cast<int32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, uint64_t &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, int64_t &a2) {
  a2 = static_cast<int64_t>(a1.val);
}

// Conversion from uint64_t

inline void TYPE_CONVERSION(stc<uint64_t> const &a1, double &a2) {
  a2 = static_cast<double>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint64_t> const &a1, uint8_t &a2) {
  a2 = static_cast<uint8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint64_t> const &a1, int8_t &a2) {
  a2 = static_cast<int8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint64_t> const &a1, uint16_t &a2) {
  a2 = static_cast<uint16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint64_t> const &a1, int16_t &a2) {
  a2 = static_cast<int16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint64_t> const &a1, uint32_t &a2) {
  a2 = static_cast<uint32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint64_t> const &a1, int32_t &a2) {
  a2 = static_cast<int32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<uint64_t> const &a1, uint64_t &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<uint64_t> const &a1, int64_t &a2) {
  a2 = static_cast<int64_t>(a1.val);
}

// size_t conversions, only when size_t is distinct from both uint64_t and
// uint32_t (e.g. Apple, where size_t is `unsigned long` while uint64_t is
// `unsigned long long`). On Linux and 32-bit, the uint64_t / uint32_t
// overloads above already cover size_t.

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<double> const &a1, T &a2) {
  a2 = a1.val;
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<int8_t> const &a1, T &a2) {
  a2 = a1.val;
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<uint8_t> const &a1, T &a2) {
  a2 = a1.val;
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<int16_t> const &a1, T &a2) {
  a2 = a1.val;
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<uint16_t> const &a1, T &a2) {
  a2 = a1.val;
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<int32_t> const &a1, T &a2) {
  a2 = a1.val;
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<uint32_t> const &a1, T &a2) {
  a2 = a1.val;
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<int64_t> const &a1, T &a2) {
  a2 = a1.val;
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1, double &a2) {
  a2 = static_cast<double>(a1.val);
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1, uint8_t &a2) {
  a2 = static_cast<uint8_t>(a1.val);
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1, int8_t &a2) {
  a2 = static_cast<int8_t>(a1.val);
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1, uint16_t &a2) {
  a2 = static_cast<uint16_t>(a1.val);
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1, int16_t &a2) {
  a2 = static_cast<int16_t>(a1.val);
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1, uint32_t &a2) {
  a2 = static_cast<uint32_t>(a1.val);
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1, int32_t &a2) {
  a2 = static_cast<int32_t>(a1.val);
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1, T &a2) {
  a2 = a1.val;
}

template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1, int64_t &a2) {
  a2 = static_cast<int64_t>(a1.val);
}

/*
  Three tiers: the identity conversion is free, an arithmetic pair without a
  dedicated TYPE_CONVERSION overload is a plain cast, and everything else
  goes through TYPE_CONVERSION as before. A dedicated overload always wins
  over the cast, so any pair with checking semantics keeps them.
 */
template <typename T1, typename T2> T1 UniversalScalarConversion(T2 const &a) {
  if constexpr (std::is_same_v<T1, T2>) {
    return a;
  } else if constexpr (!requires(stc<T2> const &x, T1 &y) {
                         TYPE_CONVERSION(x, y);
                       } && std::is_arithmetic_v<T1> &&
                       std::is_arithmetic_v<T2>) {
    return static_cast<T1>(a);
  } else {
    T1 ret;
    try {
      stc<T2> stc_a{a};
      TYPE_CONVERSION(stc_a, ret);
    } catch (ConversionException &e) {
      std::cerr << "ConversionError e=" << e.val << "\n";
      throw TerminalException{1};
    }
    return ret;
  }
}

/*
  T_frexp(x, e) returns m with x = m 2^e and 1/2 <= |m| < 1, or m = 0 and
  e = 0 for x = 0: std::frexp for any number type.

  Unlike UniversalScalarConversion<double, T> it cannot overflow. An exact
  type holds values far beyond the range of a double, and a caller needing
  only the magnitude of a quotient or a square root of such values forms it
  from the mantissas and exponents, which stay in range whenever the result
  does.

  The accuracy is that of a double for the native, integer and rational
  types, and for QuadField. For RealField and RealRing it is that of their
  evaluation at the double approximation of the generator, as for get_d: no
  overflow, but a sum that cancels loses the digits it cancels.

  The generic form below goes through a double, so it is right only for a
  type whose values fit one, which covers the native types. Every other type
  overloads it next to its double conversion. The helpers that follow
  combine numbers already in this form, so that no overload has to form a
  double of a large quantity.
 */

// ldexp with the exponent clamped to where the result is already 0 or
// infinite, so that an exponent held in a long cannot overflow the int.
inline double T_frexp_ldexp(double const &m, long const &e) {
  long const e_max = 4096;
  long e_clamp = e < -e_max ? -e_max : (e > e_max ? e_max : e);
  return std::ldexp(m, static_cast<int>(e_clamp));
}

// m 2^e for any finite double m, returned in the form of T_frexp.
inline double T_frexp_normalize(double const &m, long const &e,
                                long &e_out) {
  int k;
  double r = std::frexp(m, &k);
  e_out = (r == 0) ? 0 : e + k;
  return r;
}

// The quotient (m1 2^e1) / (m2 2^e2), m2 nonzero, in the form of T_frexp.
inline double T_frexp_quotient(double const &m1, long const &e1,
                               double const &m2, long const &e2, long &e) {
  return T_frexp_normalize(m1 / m2, e1 - e2, e);
}

// The sum m1 2^e1 + m2 2^e2 in the form of T_frexp. The mantissas need not
// be normalized. A sum that cancels loses what it cancels, as in a double.
inline double T_frexp_sum(double const &m1, long const &e1, double const &m2,
                          long const &e2, long &e) {
  if (m1 == 0) {
    return T_frexp_normalize(m2, e2, e);
  }
  if (m2 == 0) {
    return T_frexp_normalize(m1, e1, e);
  }
  long e_max = e1 > e2 ? e1 : e2;
  double s = T_frexp_ldexp(m1, e1 - e_max) + T_frexp_ldexp(m2, e2 - e_max);
  return T_frexp_normalize(s, e_max, e);
}

template <typename T> inline double T_frexp(T const &x, long &e) {
  int e_int;
  double m = std::frexp(UniversalScalarConversion<double, T>(x), &e_int);
  e = e_int;
  return m;
}

template <typename T1, typename T2>
std::optional<T1> UniversalScalarConversionCheck(T2 const &a) {
  T1 ret;
  try {
    stc<T2> stc_a{a};
    TYPE_CONVERSION(stc_a, ret);
  } catch (ConversionException &e) {
    return {};
  }
  return ret;
}

template <typename T1, typename T2>
std::vector<T1> UniversalStdVectorScalarConversion(std::vector<T2> const &V) {
  size_t len = V.size();
  std::vector<T1> V_ret(len);
  for (size_t i = 0; i < len; i++) {
    V_ret[i] = UniversalScalarConversion<T1, T2>(V[i]);
  }
  return V_ret;
}

//
// ScalingInteger that is find a positive number an integer number q =
// ScalingInteger(x) such that q x belongs to an integer ring.
// ---For x a rational this is the denominator
// ---For x in a quadratic number field, q x should belong to something like
// Z[sqrt(d)]
//

inline void ScalingInteger_Kernel([[maybe_unused]] stc<int> const &x,
                                  int &x_ret) {
  x_ret = 1;
}

inline void ScalingInteger_Kernel([[maybe_unused]] stc<long> const &x,
                                  long &x_ret) {
  x_ret = 1;
}

template <typename T1, typename T2> T1 ScalingInteger(T2 const &a) {
  T1 ret;
  stc<T2> stc_a{a};
  ScalingInteger_Kernel(stc_a, ret);
  return ret;
}

//
// Nearest / Floor / Ceil operations
//

inline void NearestInteger_double_int(double const &xI, int &xO) {
  double xRnd_d = round(xI);
  int xRnd_z = static_cast<int>(xRnd_d);
  auto GetErr = [&](int const &u) -> double {
    double diff = static_cast<double>(u) - xI;
    if (diff < 0)
      return -diff;
    return diff;
  };
  double err = GetErr(xRnd_z);
  while (true) {
    bool IsOK = true;
    for (int i = 0; i < 2; i++) {
      int shift = 2 * i - 1;
      int xTest = xRnd_z + shift;
      double TheErr = GetErr(xTest);
      if (TheErr < err) {
        IsOK = false;
        xRnd_z = xTest;
      }
    }
    if (IsOK)
      break;
  }
  xO = xRnd_z;
}

template <typename To> void NearestInteger_double_To(double const &xI, To &xO) {
  double xRnd_d = round(xI);
  int xRnd_i = static_cast<int>(xRnd_d);
  To xRnd_To = xRnd_i;
  auto GetErr = [&](To const &u) -> double {
    double u_doubl = UniversalScalarConversion<double, To>(u);
    double diff = u_doubl - xI;
    if (diff < 0)
      return -diff;
    return diff;
  };
  double err = GetErr(xRnd_To);
  while (true) {
    bool IsOK = true;
    for (int i = 0; i < 2; i++) {
      int shift = 2 * i - 1;
      To xTest = xRnd_To + shift;
      double TheErr = GetErr(xTest);
      if (TheErr < err) {
        IsOK = false;
        xRnd_To = xTest;
      }
    }
    if (IsOK)
      break;
  }
  xO = xRnd_To;
}

// clang-format off
#endif  // SRC_NUMBER_TYPECONVERSION_H_
// clang-format on
