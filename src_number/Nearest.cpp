// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheory.h"
#include "factorizations.h"
#include "Timings.h"
// clang-format on

template<typename T, typename Tint>
void process() {
  for (int i=1; i<10; i++) {
    for (int j=1; j<10; j++) {
      T x_i(i);
      T x_j(j);
      //
      T quot = x_i / x_j;
      Tint quot_n = UniversalNearestScalarInteger<Tint,T>(quot);
      std::cerr << "quot=" << quot << " quot_n" << quot_n << "\n";
    }
  }

}



int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    {
      using T = mpq_class;
      using Tint = mpz_class;
      process<T, Tint>();
    }
    {
      using T = boost::multiprecision::mpq_rational;
      using Tint = boost::multiprecision::mpz_int;
      process<T, Tint>();
    }
    {
      using T = boost::multiprecision::cpp_rational;
      using Tint = boost::multiprecision::cpp_int;
      process<T, Tint>();
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
