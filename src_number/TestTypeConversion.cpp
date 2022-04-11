#include "NumberTheoryGmp.h"

int main() {

  for (int i_int = 0; i_int < 10; i_int++) {
    uint8_t i_uint8_t = i_int;
    uint16_t i_uint16_t = i_int;
    uint32_t i_uint32_t = i_int;
    mpz_class val1 = UniversalScalarConversion<mpz_class, int>(i_int);
    mpz_class val2 = UniversalScalarConversion<mpz_class, uint8_t>(i_uint8_t);
    mpz_class val3 = UniversalScalarConversion<mpz_class, uint16_t>(i_uint16_t);
    mpz_class val4 = UniversalScalarConversion<mpz_class, uint32_t>(i_uint32_t);
    std::cerr << "val1=" << val1 << " val2=" << val2 << " val3=" << val3
              << " val4=" << val4 << "\n";
  }
}
