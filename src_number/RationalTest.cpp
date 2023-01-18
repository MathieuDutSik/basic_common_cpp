// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "rational.h"

#include "NumberTheory.h"
#include "factorizations.h"

int main(int argc, char *argv[]) {
  using T = Rational<mpz_class>;
  using Tidx_value = int;
  try {
    std::unordered_map<T, Tidx_value> ValueMap;
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
