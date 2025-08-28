// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "TypeConversion.h"
#include "TypeConversionFinal.h"
// clang-format on

// This is mostly a compilation test

template <typename T, typename Tred> void process_A() {
  T val1(23);
  Tred val2 = UniversalScalarConversion<Tred, T>(val1);
  T val3 = UniversalScalarConversion<T, Tred>(val2);
  if (val1 != val3) {
    std::cerr << "Some incoherence\n";
    throw TerminalException{1};
  }
}

template <typename T> void process_B() {
  process_A<T, int8_t>();
  process_A<T, int16_t>();
  process_A<T, int32_t>();
  process_A<T, int64_t>();
}

int main(int argc, char *argv[]) {
  try {
    process_B<mpq_class>();
    process_B<mpz_class>();
    process_B<boost::multiprecision::cpp_int>();
    process_B<boost::multiprecision::cpp_rational>();
    process_B<boost::multiprecision::mpz_int>();
    process_B<boost::multiprecision::mpq_rational>();
  } catch (TerminalException const &e) {
    std::cerr << "There is a bug to be resolved in the code\n";
    exit(e.eVal);
  }
}
