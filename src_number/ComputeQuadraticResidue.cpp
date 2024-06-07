// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include "TypeConversion.h"
#include "quadratic_residue.h"
// clang-format on

template <typename T>
void process(std::string const& a_str, std::string const& p_str, std::string const& OutFormat, std::ostream& os) {
  T a = ParseScalar<T>(a_str);
  T p = ParseScalar<T>(p_str);

  bool test = is_quadratic_residue(a, p);
  if (OutFormat == "GAP") {
    os << "return " << test << ";\n";
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "Program is used as\n";
      std::cerr << "ComputeQuadraticResidue [a] [p] [OutFormat] [FileOut]\n";
      std::cerr << "or\n";
      std::cerr << "ComputeQuadraticResidue [a] [p]\n";
      return 1;
    }
    std::string a_str = argv[1];
    std::string p_str = argv[2];
    std::string OutFormat = "GAP";
    std::string FileOut = "stderr";
    using T = mpz_class;
    if (argc == 5) {
      OutFormat = argv[3];
      FileOut = argv[4];
    }
    if (FileOut == "stderr") {
      process<T>(a_str, p_str, OutFormat, std::cerr);
    } else {
      if (FileOut == "stdout") {
        process<T>(a_str, p_str, OutFormat, std::cout);
      } else {
        std::ofstream os(FileOut);
        process<T>(a_str, p_str, OutFormat, os);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "There is a bug to be resolved in the code\n";
    exit(e.eVal);
  }
}
