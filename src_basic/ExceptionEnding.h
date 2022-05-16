// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_EXCEPTIONENDING_H_
#define SRC_BASIC_EXCEPTIONENDING_H_

#include <cassert>
#include <string>

// types for exception
struct TerminalException {
  int eVal;
};

struct ConversionException {
  std::string val;
};

// This is guaranteed to trigger an end.
// Also it gives something that can be used for having the stacktrace via gdb.
void TerminalEnding() {
  assert(false);
  exit(1);
}

// clang-format off
#endif  // SRC_BASIC_EXCEPTIONENDING_H_
// clang-format on
