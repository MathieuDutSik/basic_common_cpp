// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_EXCEPTIONSFUNC_H_
#define SRC_BASIC_EXCEPTIONSFUNC_H_

#include <cassert>
#include <cstdint>
#include <string>

// types for exception
struct TerminalException {
  int eVal;
};

struct RuntimeException {
  int eVal;
};

struct SafeIntException {
  int64_t eVal;
};

// Thrown by terminate_in_arithmetic_error<TryInt64>() when an overflow was
// recorded by the deferred-checking TryInt64 arithmetic (see
// NumberTheorySafeInt.h), and by the conversions into TryInt64 when a value
// does not fit into int64_t.
struct TryIntException {
  int64_t eVal;
};

// This is guaranteed to trigger an end.
// Also it gives something that can be used for having the stacktrace via gdb.
inline void TerminalEnding() {
  assert(false);
  exit(1);
}

// clang-format off
#endif  // SRC_BASIC_EXCEPTIONSFUNC_H_
// clang-format on
