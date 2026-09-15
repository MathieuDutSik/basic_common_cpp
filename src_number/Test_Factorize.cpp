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

template <typename T> std::string test(std::string name_numeric) {
  HumanTime time;
  std::stringstream os;
  for (int n = 1; n < 500; n++) {
    T n_T = n;
    std::vector<T> V = FactorsInt(n_T);
    // FactorsInt must return the prime factors in non-decreasing order.
    for (size_t i = 1; i < V.size(); i++) {
      if (V[i] < V[i - 1]) {
        std::cerr << "FactorsInt is not sorted for n=" << n_T << "\n";
        throw TerminalException{1};
      }
    }
    os << "n=" << n_T << " Fact=";
    for (auto &val : V)
      os << " " << val;
    os << "\n";
    std::vector<T> Ldiv = GetAllFactors(n_T);
    os << "  divisors =";
    for (auto &val : Ldiv)
      os << " " << val;
    os << "\n";
    // Also exercise FactorsIntMap (the function used by Factorize.cpp) and
    // verify it matches the multiset implied by FactorsInt.
    std::map<T, size_t> M = FactorsIntMap(n_T);
    std::map<T, size_t> M_from_V;
    for (auto &val : V)
      M_from_V[val]++;
    if (M != M_from_V) {
      std::cerr << "FactorsIntMap does not agree with FactorsInt for n="
                << n_T << "\n";
      throw TerminalException{1};
    }
    os << "  factor_map =";
    for (auto &kv : M)
      os << " " << kv.first << "^" << kv.second;
    os << "\n";
  }
  std::cerr << "Result for numeric=" << name_numeric << " time=" << time
            << "\n";
  std::string converted(os.str());
  return converted;
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2) {
      std::cerr << "Factorize [oper]\n";
      std::cerr << "\n";
      std::cerr << "oper: check the \n";
      std::cerr << "print: output the data to the file\n";
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
        std::cerr << "SafeInt64 : " << test<SafeInt64>("SafeInt64") << "\n";
        std::cerr << "cpp_int : "
                  << test<boost::multiprecision::cpp_int>("cpp_int") << "\n";
        std::cerr << "mpz_int : "
                  << test<boost::multiprecision::mpz_int>("mpz_int") << "\n";
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
