// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "TypeConversion.h"
#include "quadratic_residue.h"
// clang-format on

template <typename T> void process(int TheMod) {
  T TheMod_T(TheMod);
  int n_quad = 0;
  for (int u = 0; u < TheMod; u++) {
    T u_T(u);
    bool test = is_quadratic_residue(u_T, TheMod_T);
    if (test) {
      n_quad += 1;
    }
  }
  int p = (TheMod - 1) / 2;
  int n_quad_exact = 1 + p;
  if (n_quad != n_quad_exact) {
    std::cerr << "TheMod=" << TheMod << "\n";
    std::cerr << "n_quad=" << n_quad << " n_quad_exact=" << n_quad_exact
              << "\n";
    throw TerminalException{1};
  }
  T test_x(1);
  T test_m(2);
  bool test = is_quadratic_residue(test_x, test_m);
  if (!test) {
    std::cerr << "Failed to correctly resolve for x=1 m=2\n";
    throw TerminalException{1};
  }
}

int main() {
  try {
    int TheMod = 23;
    std::cerr << "TheMod=" << TheMod << "\n";
    process<mpz_class>(TheMod);
    process<SafeInt64>(TheMod);
    //    process<boost::multiprecision::cpp_int>(TheMod);
    //    process<boost::multiprecision::mpz_int>(TheMod);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "There is a bug to be resolved in the code\n";
    exit(e.eVal);
  }
}
