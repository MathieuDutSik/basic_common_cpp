// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "rational.h"
#include "NumberTheoryCommon.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "factorizations.h"
// clang-format on

int main() {
  using Tidx_value = int;
  try {
    std::unordered_map<Rational<mpz_class>, Tidx_value> ValueMap1;
    std::unordered_map<Rational<SafeInt64>, Tidx_value> ValueMap2;
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
