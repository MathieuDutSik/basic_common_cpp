// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "factorizations.h"
#include "SafeInteger.h"
// clang-format on

template<typename T>
void test() {
  for (int n = 1; n < 500; n++) {
    T n_T = n;
    std::vector<T> V = FactorsInt(n_T);
    std::cerr << "n=" << n_T << " Fact=";
    for (auto &val : V)
      std::cerr << " " << val;
    std::cerr << "\n";
    std::vector<T> Ldiv = GetAllFactors(n_T);
    std::cerr << "  divisors =";
    for (auto &val : Ldiv)
      std::cerr << " " << val;
    std::cerr << "\n";
  }
}



int main() {
  try {
    test<mpz_class>();
    test<mpq_class>();
    test<SafeInt64>();
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
