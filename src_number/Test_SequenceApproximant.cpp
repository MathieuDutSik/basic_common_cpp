// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheory.h"
#include "fractions.h"
#include "Timings.h"
#include "NumberTheorySafeInt.h"
// clang-format on

template <typename T> std::string test(std::string name_numeric) {
  HumanTime time;
  std::stringstream os;
  for (int n = 1; n < 100; n++) {
    int num = rand() % 10000;
    int den = rand() % 10000;
    int choice = rand() % 2;
    int sign = 2 * choice - 1;
    T num_T(num);
    T den_T(den);
    T val = sign * num_T / den_T;
    std::vector<T> list_approx =
        get_sequence_continuous_fraction_approximant(val);
    os << "val=" << val << " list_approx =";
    for (auto &approx : list_approx) {
      os << " " << approx;
    }
    os << "\n";
    T last_approx = list_approx[list_approx.size() - 1];
    if (last_approx != val) {
      std::cerr << "last_approx=" << last_approx << " val=" << val << "\n";
      throw TerminalException{1};
    }
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
      std::cerr << "test_sequence_approximant [oper]\n";
      std::cerr << "\n";
      std::cerr << "oper: check the \n";
      std::cerr << "print: output the data to the file\n";
      throw TerminalException{1};
    }
    std::string oper = argv[1];
    auto eval = [&]() -> void {
      if (oper == "check") {
        std::unordered_map<std::string, std::string> map;
        map[test<mpq_class>("mpz_class")] = "mpq_class";
        //        map[test<boost::multiprecision::mpq_rational>("mpq_rationql")]
        //        = "boost::mpq_rational";
        //        map[test<boost::multiprecision::cpp_rational>("cpp_rational")]
        //        = "boost::cpp_rational";
        map[test<Rational<SafeInt64>>("Rational<SafeInt64>")] =
            "Rational<SafeInt64>";
        if (map.size() != 1) {
          std::cerr << "We have incoherent result for arithmetics\n";
          std::cerr << "|map|=" << map.size() << "\n";
          throw TerminalException{1};
        }
        return;
      }
      if (oper == "print") {
        std::cerr << "mpq_class : " << test<mpq_class>("mpq_class") << "\n";
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
