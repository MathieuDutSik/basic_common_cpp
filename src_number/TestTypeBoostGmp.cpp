// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberTheoryBoostGmpInt.h"

int main() {

  using T = boost::multiprecision::mpq_rational;
  using Tint = boost::multiprecision::mpz_int;

  T val1, val2;
  val1 = 2;
  val2 = 5;
  T val = val1 / val2;
  std::cerr << "val=" << val << "\n";
  Tint val_red = UniversalNearestScalarInteger<Tint,T>(val);
  std::cerr << "val_red=" << val_red << "\n";
  //
  Tint val3 = 15;
  Tint val4 = 12;
  Tint lcm = LCMpair(val3, val4);
  std::cerr << "lcm=" << lcm << "\n";
}
