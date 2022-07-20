// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheoryGmp.h"

int main() {

  for (int i_int = 0; i_int < 10; i_int++) {
    uint8_t i_uint8_t = i_int;
    uint16_t i_uint16_t = i_int;
    uint32_t i_uint32_t = i_int;
    size_t val_sz = i_int;
    mpz_class val1 = UniversalScalarConversion<mpz_class, int>(i_int);
    mpz_class val2 = UniversalScalarConversion<mpz_class, uint8_t>(i_uint8_t);
    mpz_class val3 = UniversalScalarConversion<mpz_class, uint16_t>(i_uint16_t);
    mpz_class val4 = UniversalScalarConversion<mpz_class, uint32_t>(i_uint32_t);
    //    uint64_t val_64 = UniversalScalarConversion<uint64_t, size_t>(val_sz);
    if (std::is_same_v<size_t,uint64_t>) {
      std::cerr << "size_t and uint64_t have the same type\n";
    } else {
      std::cerr << "size_t and uint64_t DO NOT have the same type\n";
    }
    std::cerr << "val1=" << val1 << " val2=" << val2 << " val3=" << val3
              << " val4=" << val4 << "\n";
  }
}
