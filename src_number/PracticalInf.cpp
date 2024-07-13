// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "factorizations.h"
// clang-format on

int main() {
  try {
    std::cerr << "Inf(mpz_class)=" << practical_infinity<mpz_class>() << "\n";
    std::cerr << "Inf(mpq_class)=" << practical_infinity<mpq_class>() << "\n";
    std::cerr << "Inf(size_t)=" << practical_infinity<size_t>() << "\n";
    std::cerr << "Inf(int)=" << practical_infinity<int>() << "\n";
    std::cerr << "Inf(uint8_t)=" << size_t(practical_infinity<uint8_t>())
              << "\n";
    std::cerr << "Inf(int8_t)=" << size_t(practical_infinity<int8_t>()) << "\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
