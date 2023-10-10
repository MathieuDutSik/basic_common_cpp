// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "TypeConversion.h"
// clang-format on

int main(int argc, char *argv[]) {
  //  using T=int;
  //  using T=long;
  //  using T=long;
  using T = boost::multiprecision::cpp_int;
  //  using T=mpq_class;
  //  using Tint=mpz_class;
  try {
    T val = 63;
    T val2 = sqrt(val);
    std::cerr << "val=" << val << " val2=" << val2 << "\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
