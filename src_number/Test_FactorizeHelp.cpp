// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "factorizations.h"
#include "Timings.h"
// clang-format on

template <typename T> std::string test(std::string const &name_numeric) {
  HumanTime time;
  std::stringstream os;
  // For each composite N, supply a partial list of helper factors and
  // verify that FactorsIntMap_help returns the same map as FactorsIntMap.
  for (int n = 2; n < 200; n++) {
    T n_T(n);
    std::map<T, size_t> ref = FactorsIntMap(n_T);

    // Three variants of the help list:
    //  - empty (should still produce the correct factorization)
    //  - the smallest prime factor of n
    //  - the smallest prime factor squared (a guaranteed divisor when the
    //    multiplicity is at least 2; otherwise harmless)
    std::vector<std::vector<T>> help_variants;
    help_variants.push_back({});
    if (!ref.empty()) {
      T smallest = ref.begin()->first;
      help_variants.push_back({smallest});
      help_variants.push_back({smallest, smallest * smallest});
    }
    for (size_t v = 0; v < help_variants.size(); v++) {
      std::map<T, size_t> got = FactorsIntMap_help(n_T, help_variants[v]);
      if (got != ref) {
        std::cerr << "FactorsIntMap_help disagrees with FactorsIntMap for n="
                  << n_T << " variant=" << v << "\n";
        throw TerminalException{1};
      }
      os << "n=" << n_T << " help_variant=" << v << " factors=";
      for (auto &kv : got) {
        os << " " << kv.first << "^" << kv.second;
      }
      os << "\n";
    }
  }
  std::cerr << "Result for numeric=" << name_numeric << " time=" << time
            << "\n";
  return os.str();
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2) {
      std::cerr << "Test_FactorizeHelp [oper]\n";
      std::cerr << "\n";
      std::cerr << "oper: check / print\n";
      throw TerminalException{1};
    }
    std::string oper = argv[1];
    auto eval = [&]() -> void {
      if (oper == "check") {
        std::unordered_map<std::string, std::string> map;
        map[test<mpz_class>("mpz_class")] = "mpz_class";
        map[test<SafeInt64>("SafeInt64")] = "SafeInt64";
        map[test<boost::multiprecision::cpp_int>("cpp_int")] = "cpp_int";
        map[test<boost::multiprecision::mpz_int>("mpz_int")] = "mpz_int";
        if (map.size() != 1) {
          std::cerr << "We have incoherent result for arithmetics\n";
          std::cerr << "|map|=" << map.size() << "\n";
          throw TerminalException{1};
        }
        return;
      }
      if (oper == "print") {
        std::cerr << "mpz_class : " << test<mpz_class>("mpz_class") << "\n";
        return;
      }
      std::cerr << "Failed to find matching oper\n";
      throw TerminalException{1};
    };
    eval();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
