// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_QUOINTFCTS_H_
#define SRC_NUMBER_QUOINTFCTS_H_

// clang-format off
#include "ExceptionsFunc.h"
#include <iostream>
// clang-format on

template <typename T> T QuoInt(T const &a, T const &b) {
  T ret;
  try {
    stc<T> stc_a{a};
    stc<T> stc_b{b};
    QUO_INT(stc_a, stc_b, ret);
  } catch (QuoIntException &e) {
    std::cerr << "QuoIntError e=" << e.val << "\n";
    throw TerminalException{1};
  }
  return ret;
}

// clang-format off
#endif  // SRC_NUMBER_QUOINTFCTS_H_
// clang-format on
