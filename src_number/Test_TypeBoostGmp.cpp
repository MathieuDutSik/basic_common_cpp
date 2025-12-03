// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "ResidueQuotient.h"
// clang-format on

template<typename Trat, typename Tint>
void test_trat_tint() {
  Trat val1, val2;
  val1 = 2;
  val2 = 5;
  Trat val = val1 / val2;
  std::cerr << "val=" << val << "\n";
  Tint val_red = UniversalNearestScalarInteger<Tint, Trat>(val);
  std::cerr << "val_red=" << val_red << "\n";
  //
  using Tlong = int64_t;
  Tlong val3 = 15;
  Tlong val4 = 12;
  Tlong lcm = LCMpair(val3, val4);
  std::cerr << "lcm=" << lcm << "\n";
}



int main() {
  {
    using Trat = boost::multiprecision::mpq_rational;
    using Tint = boost::multiprecision::mpz_int;
    test_trat_tint<Trat,Tint>();
  }
  //
  {
    using Trat = boost::multiprecision::cpp_rational;
    using Tint = boost::multiprecision::cpp_int;
    test_trat_tint<Trat,Tint>();
  }
}
