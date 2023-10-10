// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_TYPECONVERSION_H_
#define SRC_NUMBER_TYPECONVERSION_H_

// clang-format off
#include "BasicNumberTypes.h"
#include "ExceptionsFunc.h"
#include "TemplateTraits.h"
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

inline void TYPE_CONVERSION(stc<double> const &a1, T_uint64_t &a2) {
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

inline void TYPE_CONVERSION(stc<int8_t> const &a1, T_uint64_t &a2) {
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

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, T_uint64_t &a2) {
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

inline void TYPE_CONVERSION(stc<int16_t> const &a1, T_uint64_t &a2) {
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

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, T_uint64_t &a2) {
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

inline void TYPE_CONVERSION(stc<int32_t> const &a1, T_uint64_t &a2) {
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

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, T_uint64_t &a2) {
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

inline void TYPE_CONVERSION(stc<int64_t> const &a1, T_uint64_t &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, int64_t &a2) {
  a2 = static_cast<int64_t>(a1.val);
}

// Conversion from uint64_t

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, double &a2) {
  a2 = static_cast<double>(a1.val);
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, uint8_t &a2) {
  a2 = static_cast<uint8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, int8_t &a2) {
  a2 = static_cast<int8_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, uint16_t &a2) {
  a2 = static_cast<uint16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, int16_t &a2) {
  a2 = static_cast<int16_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, uint32_t &a2) {
  a2 = static_cast<uint32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, int32_t &a2) {
  a2 = static_cast<int32_t>(a1.val);
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, T_uint64_t &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, int64_t &a2) {
  a2 = static_cast<int64_t>(a1.val);
}

template <typename T1, typename T2> T1 UniversalScalarConversion(T2 const &a) {
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

void NearestInteger_double_int(double const &xI, int &xO) {
  //  std::cerr << "Temp_common : NearestInteger\n";
  double xRnd_d = round(xI);
  int xRnd_z = static_cast<int>(xRnd_d);
  //  std::cerr << "xI=" << xI << "\n";
  auto GetErr = [&](int const &u) -> double {
    double diff = static_cast<double>(u) - xI;
    if (diff < 0)
      return -diff;
    return diff;
  };
  double err = GetErr(xRnd_z);
  //  std::cerr << "err=" << err << "\n";
  while (true) {
    bool IsOK = true;
    for (int i = 0; i < 2; i++) {
      int shift = 2 * i - 1;
      int xTest = xRnd_z + shift;
      double TheErr = GetErr(xTest);
      //      std::cerr << "i=" << i << " shift=" << shift << " xTest=" << xTest
      //      << " TheErr=" << TheErr << "\n";
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
  //  std::cerr << "Temp_common : NearestInteger\n";
  double xRnd_d = round(xI);
  int xRnd_i = static_cast<int>(xRnd_d);
  To xRnd_To = xRnd_i;
  //  std::cerr << "xI=" << xI << "\n";
  auto GetErr = [&](To const &u) -> double {
    double u_doubl = UniversalScalarConversion<double, To>(u);
    double diff = u_doubl - xI;
    if (diff < 0)
      return -diff;
    return diff;
  };
  double err = GetErr(xRnd_To);
  //  std::cerr << "err=" << err << "\n";
  while (true) {
    bool IsOK = true;
    for (int i = 0; i < 2; i++) {
      int shift = 2 * i - 1;
      To xTest = xRnd_To + shift;
      double TheErr = GetErr(xTest);
      //      std::cerr << "i=" << i << " shift=" << shift << " xTest=" << xTest
      //      << " TheErr=" << TheErr << "\n";
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
