// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheory.h"
#include "factorizations.h"
// clang-format on

template <typename Trat> void process(std::string const &name) {
  using T = QuadField<Trat, 5>;
  T x;
  T y = UniversalScalarConversion<T, T>(x);
  std::cerr << "name=" << name << " x=" << x << " y=" << y << "\n";
  Trat half = Trat(1) / Trat(2);
  T phi(half, half);
  T pow(1);
  for (int i = 0; i < 10; i++) {
    std::cerr << "i=" << i << " pow=" << pow << "\n";
    pow *= phi;
  }
  T near = UniversalNearestScalarInteger<T, T>(x);
  std::cerr << "near=" << near << "\n";
}

int main() {
  try {
    process<mpq_class>("mpq_class");
    process<Rational<SafeInt64>>("Rational<SafeInt64>");
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
