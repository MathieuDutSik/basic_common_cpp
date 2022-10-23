// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "MAT_MatrixInt.h"
#include "NumberTheory.h"
#include "QuadField.h"
int main(int argc, char *argv[]) {
  using Trat = mpq_class;
  using T = QuadField<Trat,5>;
  try {
    MyVector<T> V(5);
    FractionVector<T> frv = NonUniqueScaleToIntegerVectorPlusCoeff(V);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
