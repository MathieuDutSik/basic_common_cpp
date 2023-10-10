// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "rational.h"
#include "NumberTheory.h"
#include "factorizations.h"
// clang-format on

int main(int argc, char *argv[]) {
  using T = Rational<mpz_class>;
  using Tidx_value = int;
  try {
    std::unordered_map<T, Tidx_value> ValueMap;
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
