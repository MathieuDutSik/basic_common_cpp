// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheoryRealField.h"
#include "NumberTheory.h"


// The quantity is 2*cos(2*pi/7)
// The minimal polynomial is X^3 + X^2 - 2X - 1

int main(int argc, char *argv[]) {
  try {
    using T_rat = mpq_class;
    std::string eFile = "../Examples/CubicField/CubicFieldDisc_49";
    HelperClassRealField<T_rat> hcrf(eFile);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
