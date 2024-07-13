// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_BASICNUMBERTYPES_H_
#define SRC_NUMBER_BASICNUMBERTYPES_H_
// clang-format off
#include <string>
#include <cstdint>
// clang-format on

// STC: Singleton Type Conversion
// We absolutely want to avoid a function with a signature "long"
// matching an int. That is why we introduce the stc<T> data type
// since C++ will never convert a stc<long> to a stc<int> and
// vice versa under the hood.
// The overhead is eliminated at the compilation. The stc<T> does
// not show up outside of internal conversion code.
template <typename T> struct stc {
  T const &val;
};

struct ConversionException {
  std::string val;
};

struct QuoIntException {
  std::string val;
};

template <typename T> struct PairGCD_dot {
  T a;
  T b;
  T gcd;
};

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

// clang-format off
#endif  // SRC_NUMBER_BASICNUMBERTYPES_H_
// clang-format on
