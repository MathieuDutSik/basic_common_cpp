// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryQuadField.h"
#include "MAT_MatrixInt.h"
// clang-format on

int main() {
  using Trat = mpq_class;
  using T = QuadField<Trat, 5>;
  try {
    MyVector<T> V(5);
    FractionVector<T> frv = NonUniqueScaleToIntegerVectorPlusCoeff(V);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
