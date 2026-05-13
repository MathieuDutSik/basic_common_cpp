// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheory.h"
// clang-format on


int main() {
  using T1 = mpq_class;
  using T1_int = mpz_class;
  using T2 = boost::multiprecision::mpq_rational;
  using T2_int = boost::multiprecision::mpz_int;
  using T3 = boost::multiprecision::cpp_rational;
  using T3_int = boost::multiprecision::cpp_int;


  size_t val = 3;
  T3_int val_b = UniversalScalarConversion<T3_int,size_t>(val);
}
