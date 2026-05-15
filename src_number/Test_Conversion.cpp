// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "TypeConversion.h"
#include "TypeConversionFinal.h"
// clang-format on

int main([[maybe_unused]] int argc, [[maybe_unused]] char *argv[]) {
  double val1_d = 1.0;
  mpz_class val1_z = UniversalScalarConversion<mpz_class, double>(val1_d);
  SafeInt64 val1_si = UniversalScalarConversion<SafeInt64, double>(val1_d);
  std::cerr << "val1_z=" << val1_z << " val1_si=" << val1_si
            << " val1_d=" << val1_d << "\n";
}
