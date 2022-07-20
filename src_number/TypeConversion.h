// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_TYPECONVERSION_H_
#define SRC_NUMBER_TYPECONVERSION_H_

// All the definitions of special fields are in other include.
// Nothing of this should depend on GMP or MPREAL or FLINT or whatever.
//
//  All GMP are in NumberTheory.h
//  All mpreal are in mpreal_related.h

#include "ExceptionEnding.h"
#include "TemplateTraits.h"
#include <cstdint>
#include <iostream>
#include <math.h>
#include <type_traits>
#include <utility>

//
// UniversalScalarConversion and TYPE_CONVERSION
//

// STC: Singleton Type Conversion
// We absolutely want to avoid a function with a signature "long"
// matching an int. That is why we introduce the stc<T> data type
// since C++ will never convert a stc<long> to a stc<int> and
// vice versa under the hood.
// The overhead is eliminated at the compilation. The stc<T> does
// not show up outside of internal conversion code.
template <typename T> struct stc { T const &val; };


// The problem we face is that we have
// --- std::is_same_v<size_t,uint64_t> = T on the Linux X86
// --- std::is_same_v<size_t,uint64_t> = F on the Macintosh X86
// By introducing that type we guarantee that we have always
// std::is_same_v<size_t,T_uint64_t> = T on both platforms
#ifdef __APPLE__
using T_uint64_t = unsigned long;
#else
using T_uint64_t = uint64_t;
#endif

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
std::pair<bool, T1> UniversalScalarConversionCheck(T2 const &a) {
  T1 ret;
  try {
    stc<T2> stc_a{a};
    TYPE_CONVERSION(stc_a, ret);
  } catch (ConversionException &e) {
    return {false, ret};
  }
  return {true, ret};
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

template <typename To, typename Ti>
inline To UniversalFloorScalarInteger(Ti const &a) {
  To ret;
  FloorInteger(a, ret);
  return ret;
}

template <typename To, typename Ti>
inline To UniversalCeilScalarInteger(Ti const &a) {
  To ret;
  CeilInteger(a, ret);
  return ret;
}

template <typename To, typename Ti>
inline typename std::enable_if<!std::is_same_v<Ti, double>, To>::type
UniversalNearestScalarInteger(Ti const &a) {
  To ret;
  NearestInteger(a, ret);
  return ret;
}

template <typename To, typename Ti>
inline typename std::enable_if<std::is_same_v<Ti, double>, To>::type
UniversalNearestScalarInteger(Ti const &a) {
  To ret;
  NearestInteger_double_To<To>(a, ret);
  return ret;
}

// clang-format off
#endif  // SRC_NUMBER_TYPECONVERSION_H_
// clang-format on
