// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_BASICNUMBERTYPES_H_
#define SRC_NUMBER_BASICNUMBERTYPES_H_
// clang-format off
#include <format>
#include <sstream>
#include <string>
#include <cstdint>
// clang-format on

/*
  Rendering of the number types through std::format.

  A specialization of std::formatter is the way the standard offers for making
  a type printable, and it is allowed: it is a specialization of a standard
  class template that depends on a program-defined type, which is what
  [namespace.std] permits. Adding overloads of std::to_string is not allowed by
  anything, whatever the argument type, since the permission covers class
  template specializations and not functions.

  Every number type of this project already has an operator<<, so the
  formatters all delegate to it. A type gains std::format support with

    template <> struct std::formatter<MyType> : ostream_formatter<MyType> {};

  and the callers use std::format("{}", x) instead of std::to_string(x).
 */
template <typename T> struct ostream_formatter {
  constexpr auto parse(std::format_parse_context &ctx) const {
    // No format specification is accepted beyond the empty one.
    return ctx.begin();
  }
  template <typename FormatContext>
  auto format(T const &val, FormatContext &ctx) const {
    std::ostringstream os;
    os << val;
    return std::format_to(ctx.out(), "{}", os.str());
  }
};

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

// clang-format off
#endif  // SRC_NUMBER_BASICNUMBERTYPES_H_
// clang-format on
