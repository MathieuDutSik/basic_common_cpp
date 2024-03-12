// Copyright (C) 2024 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_FRACTIONS_H_
#define SRC_NUMBER_FRACTIONS_H_

// clang-format off
#include "TemplateTraits.h"
#include <vector>
// clang-format on

template <typename T>
std::vector<T> get_continous_fraction(T const& x) {
  static_assert(is_implementation_of_Q<T>::value, "Requires T to be a rational field");
  std::vector<T> approx;
  T x_work = x;
  while (true) {
    T low_val = UniversalFloorScalarInteger<T,T>(x_work);
    approx.push_back(low_val);
    if (low_val == x_work) {
      return approx;
    }
    x_work = 1 / (x_work - low_val);
  }
}

template<typename T>
std::vector<T> get_sequence_continuous_fraction_approximant(T const& x) {
  std::vector<T> approx = get_continous_fraction(x);
  size_t len = approx.size();
  std::vector<T> list_approx;
  for (size_t u=0; u<len; u++) {
    T sing_approx = approx[u];
    if (u > 0) {
      for (size_t pos1=0; pos1<=u-1; pos1++) {
        size_t pos = u - 1 - pos1;
        T val = approx[pos];
        sing_approx = val + 1/sing_approx;
      }
    }
    list_approx.push_back(sing_approx);
  }
  return list_approx;
}

template<typename T>
T get_mid_val(T const& TheLow, T const& TheUpp) {
  T eFrac = (TheLow + TheUpp) / 2;
  T TargetLow = (2*TheLow + TheUpp) / 3;
  T TargetUpp = (TheLow + 2*TheUpp) / 3;
  std::vector<T> list_approx = get_sequence_continuous_fraction_approximant(eFrac);
  for (auto & approx : list_approx) {
    if (TargetLow <= approx && approx <= TargetUpp) {
      return approx;
    }
  }
  std::cerr << "We should not reach that state in GetMidVal\n";
  throw TerminalException{1};
};

// clang-format off
#endif  // SRC_NUMBER_FRACTIONS_H_
// clang-format on
