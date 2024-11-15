// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "factorizations.h"
#include "Timings.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2) {
      std::cerr << "Factorize [N]\n";
      std::cerr << "\n";
      std::cerr << "N: check the \n";
      std::cerr << "print: output the data to the file\n";
      throw TerminalException{1};
    }
    using T = mpz_class;
    std::string N_str = argv[1];
    T N = ParseScalar<T>(N_str);
    std::map<T, size_t> map = FactorsIntMap(N);
    for (auto & kv : map) {
      std::cerr << "p=" << kv.first << " mult=" << kv.second << "\n";
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
