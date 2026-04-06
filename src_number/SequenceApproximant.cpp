// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheory.h"
#include "fractions.h"
#include "Timings.h"
#include "NumberTheorySafeInt.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2) {
      std::cerr << "SequenceApproximant [number]\n";
      throw TerminalException{1};
    }
    std::string input_s = argv[1];
    using T = mpq_class;
    T input = ParseScalar<T>(input_s);

    std::vector<T> list_approx =
        get_sequence_continuous_fraction_approximant(input);
    size_t len = list_approx.size();
    for (size_t i=0; i<len; i++) {
      T approx = list_approx[i];
      double approx_d = UniversalScalarConversion<double,T>(approx);
      std::cerr << "i=" << i << "/" << len << " approx=" << approx << " approx_d=" << approx_d << "\n";
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
