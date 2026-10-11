// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheory.h"
// clang-format on


int main() {
  using T3_int = boost::multiprecision::cpp_int;
  size_t val = 3;
  T3_int val_b = UniversalScalarConversion<T3_int, size_t>(val);
  std::cerr << "val=" << val << " val_b=" << val_b << "\n";
}
