// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryGmp.h"
// clang-format on

int main() {
  std::cerr << "sizeof(short int)=" << sizeof(short int) << "\n";
  std::cerr << "sizeof(int)=" << sizeof(int) << "\n";
  std::cerr << "sizeof(int8_t)=" << sizeof(int8_t) << "\n";
  std::cerr << "sizeof(int64_t)=" << sizeof(int64_t) << "\n";
  std::cerr << "sizeof(long)=" << sizeof(long) << "\n";
  std::cerr << "sizeof(long long)=" << sizeof(long long) << "\n";
  if (std::is_same_v<short int, int16_t>) {
    std::cerr << "short int and int16_t have the same type\n";
  } else {
    std::cerr << "short int and int16_t DO NOT have the same type\n";
  }
  if (std::is_same_v<int, int32_t>) {
    std::cerr << "int and int32_t have the same type\n";
  } else {
    std::cerr << "int and int32_t DO NOT have the same type\n";
  }
  if (std::is_same_v<long, int64_t>) {
    std::cerr << "long and int64_t have the same type\n";
  } else {
    std::cerr << "long and int64_t DO NOT have the same type\n";
  }
  //
  if (std::is_same_v<size_t, uint64_t>) {
    std::cerr << "size_t and uint64_t have the same type\n";
  } else {
    std::cerr << "size_t and uint64_t DO NOT have the same type\n";
  }
  //
  if (std::is_same_v<size_t, T_uint64_t>) {
    std::cerr << "size_t and T_uint64_t have the same type\n";
  } else {
    std::cerr << "size_t and T_uint64_t DO NOT have the same type\n";
  }
  //
  if (std::is_same_v<size_t, unsigned long>) {
    std::cerr << "size_t and unsigned long have the same type\n";
  } else {
    std::cerr << "size_t and unsigned long DO NOT have the same type\n";
  }
}
