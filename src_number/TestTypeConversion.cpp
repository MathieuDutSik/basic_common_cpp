// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheoryGmp.h"

int main() {

  for (int i_int = 0; i_int < 10; i_int++) {
    uint8_t i_uint8_t = i_int;
    uint16_t i_uint16_t = i_int;
    uint32_t i_uint32_t = i_int;
    T_uint64_t val_uint64_t = i_int;
    size_t val_sz = i_int;
    mpz_class val_mpz = i_int;
    mpz_class val1 = UniversalScalarConversion<mpz_class, int>(i_int);
    mpz_class val2 = UniversalScalarConversion<mpz_class, uint8_t>(i_uint8_t);
    mpz_class val3 = UniversalScalarConversion<mpz_class, uint16_t>(i_uint16_t);
    mpz_class val4 = UniversalScalarConversion<mpz_class, uint32_t>(i_uint32_t);
    T_uint64_t val_64_1 = UniversalScalarConversion<T_uint64_t, size_t>(val_sz);
    T_uint64_t val_64_2 =
        UniversalScalarConversion<T_uint64_t, mpz_class>(val_mpz);
    size_t val_64_3 = UniversalScalarConversion<size_t, mpz_class>(val_mpz);
    mpz_class val_mpz_1 = UniversalScalarConversion<mpz_class, size_t>(val_sz);
    mpz_class val_mpz_2 =
        UniversalScalarConversion<mpz_class, T_uint64_t>(val_uint64_t);

    if (std::is_same_v<size_t, uint64_t>) {
      std::cerr << "size_t and uint64_t have the same type\n";
    } else {
      std::cerr << "size_t and uint64_t DO NOT have the same type\n";
    }
    if (std::is_same_v<size_t, T_uint64_t>) {
      std::cerr << "size_t and T_uint64_t have the same type\n";
    } else {
      std::cerr << "size_t and T_uint64_t DO NOT have the same type\n";
    }
    if (std::is_same_v<size_t, unsigned long>) {
      std::cerr << "size_t and unsigned long have the same type\n";
    } else {
      std::cerr << "size_t and unsigned long DO NOT have the same type\n";
    }
    std::cerr << "val1=" << val1 << " val2=" << val2 << " val3=" << val3
              << " val4=" << val4 << "\n";
  }
}
