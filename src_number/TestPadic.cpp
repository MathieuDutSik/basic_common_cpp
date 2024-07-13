// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "TypeConversion.h"
// clang-format on

int main(int argc, char *argv[]) {
  using T = mpq_class;
  try {
    T val = -12;

  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
