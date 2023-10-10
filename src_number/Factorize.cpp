// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "factorizations.h"
#include "Timings.h"
#include "NumberTheorySafeInt.h"
// clang-format on

template <typename T> std::string test(std::string name_numeric) {
  HumanTime time;
  std::stringstream os;
  for (int n = 1; n < 500; n++) {
    T n_T = n;
    std::vector<T> V = FactorsInt(n_T);
    os << "n=" << n_T << " Fact=";
    for (auto &val : V)
      os << " " << val;
    os << "\n";
    std::vector<T> Ldiv = GetAllFactors(n_T);
    os << "  divisors =";
    for (auto &val : Ldiv)
      os << " " << val;
    os << "\n";
  }
  std::cerr << "Result for numeric=" << name_numeric << " time=" << time
            << "\n";
  std::string converted(os.str());
  return converted;
}

int main(int argc, char *argv[]) {
  SingletonTime time1;
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
        map[test<mpq_class>("mpq_class")] = "mpq_class";
        map[test<SafeInt64>("SafeInt64")] = "SafeInt64";
        if (map.size() != 1) {
          std::cerr << "We have incoherent result for arithmetics\n";
          std::cerr << "|map|=" << map.size() << "\n";
          throw TerminalException{1};
        }
        return;
      }
      if (oper == "print") {
        std::cerr << "mpz_class : " << test<mpz_class>("mpz_class") << "\n";
        std::cerr << "mpq_class : " << test<mpq_class>("mpq_class") << "\n";
        std::cerr << "SafeInt64 : " << test<SafeInt64>("SafeInt64") << "\n";
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
  runtime(time1);
}
