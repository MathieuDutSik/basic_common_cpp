// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "TypeConversion.h"
// clang-format on

template <typename T> void process(std::string const &name) {
  T val(64);
  std::optional<T> val2 = UniversalSquareRoot(val);
  if (val2) {
    std::cerr << "name=" << name << " val=" << val << " sqrt=" << *val2
              << "\n";
  } else {
    std::cerr << "name=" << name << " val=" << val << " is not a square\n";
  }
}

int main() {
  try {
    process<boost::multiprecision::cpp_int>("cpp_int");
    process<mpz_class>("mpz_class");
    process<SafeInt64>("SafeInt64");
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
