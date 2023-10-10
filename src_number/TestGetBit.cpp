// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
// clang-format on

template <typename T> void process_single(T const &val_t) {
  SafeInt64 val_64(val_t);
  mpz_class val_z = val_t;
  std::cerr << "mpz_class(val_t)=" << val_z << "\n";
  std::cerr << "get_bit(val_z)=" << get_bit(val_z)
            << " get_bit(-val_z)=" << get_bit(-val_z) << "\n";
  size_t siz_t = get_bit(val_t);
  std::cerr << "We have siz_t\n";
  size_t siz_64 = get_bit(val_64);
  std::cerr << "We have siz_64\n";
  size_t siz_z = get_bit(val_z);
  std::cerr << "We have siz_z\n";
  if (siz_z != siz_64 || siz_t != siz_z) {
    std::cerr << "val_t=" << val_t << " val_64=" << val_64 << " val_z=" << val_z
              << "\n";
    std::cerr << "siz_t=" << siz_t << " siz_64=" << siz_64 << " siz_z=" << siz_z
              << "\n";
    throw TerminalException{1};
  }
}

template <typename T> void process_all() {
  // The extreme value creates us trouble, so let's avoid them.
  T val_min = std::numeric_limits<T>::min() + 1;
  T val_max = std::numeric_limits<T>::max() - 1;
  for (T x = val_min; x <= val_max; x++)
    process_single<T>(x);
}

int main(int argc, char *argv[]) {
  try {
    process_all<int8_t>();
    process_all<int16_t>();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
}
