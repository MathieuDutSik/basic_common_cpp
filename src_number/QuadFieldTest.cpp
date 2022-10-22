// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "QuadField.h"
#include "NumberTheory.h"
#include "factorizations.h"


int main(int argc, char *argv[]) {
  using T = QuadField<mpq_class,5>;
  try {
    T x;
    T y = UniversalScalarConversion<T,T>(x);
    std::cerr << "x=" << x << " y=" << y << "\n";

  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
