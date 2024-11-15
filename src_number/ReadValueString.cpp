// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "InputOutput.h"
#include "Timings.h"
// clang-format on

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2) {
      std::cerr << "ReadValueString [value]\n";
      throw TerminalException{1};
    }
    using T = mpq_class;
    std::string input = argv[1];
    size_t deg = 4;
    std::istringstream is(input);
    std::vector<T> v = ReadVectorFromRealAlgebraicString<T>(is, deg);
    for (size_t u = 0; u < deg; u++) {
      std::cerr << "u=" << u << " v[u]=" << v[u] << "\n";
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
