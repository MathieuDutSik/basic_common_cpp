// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
// clang-format on

int main(int argc, char *argv[]) {
  std::cerr << "sizeof(int8_t)=" << sizeof(int8_t) << "\n";
  std::cerr << "sizeof(int64_t)=" << sizeof(int64_t) << "\n";
  std::cerr << "sizeof(long)=" << sizeof(long) << "\n";
  std::cerr << "sizeof(long long)=" << sizeof(long long) << "\n";
}
