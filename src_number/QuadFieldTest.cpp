// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "QuadField.h"
#include "NumberTheory.h"
#include "factorizations.h"

int main(int argc, char *argv[]) {
  using Trat = mpq_class;
  using T = QuadField<Trat, 5>;
  try {
    T x;
    T y = UniversalScalarConversion<T, T>(x);
    std::cerr << "x=" << x << " y=" << y << "\n";
    Trat half = Trat(1) / Trat(2);
    T phi(half, half);
    T pow(1);
    for (int i = 0; i < 10; i++) {
      std::cerr << "i=" << i << " pow=" << pow << "\n";
      pow *= phi;
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
