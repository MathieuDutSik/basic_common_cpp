// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_APPROXIMATION_H_
#define SRC_NUMBER_APPROXIMATION_H_

#include <map>
#include <utility>
#include <vector>

// This code is used for checking the numerics comparison code
// It is not really useful since continued fractions work better
template <typename T>
T find_approximation_dichotomy(T const& val, T const& thr) {
  T val_low, val_upp;
  int pos = 1;
  while(true) {
    val_low = -pos;
    val_upp = pos;
    if (val_low < val && val < val_upp)
      break;
    pos++;
  }
  // We now have our starting point now reducing it.
  T mid;
  while(true) {
    mid = (val_upp + val_low) / 2;
    T delta = val_upp - val_low;
    if (delta < thr)
      break;
    if (val < mid) {
      val_upp = mid;
    } else {
      if (val > mid)
        val_low = mid;
    }
  }
  return mid;
}

// clang-format off
#endif  // SRC_NUMBER_APPROXIMATION_H_
// clang-format on
