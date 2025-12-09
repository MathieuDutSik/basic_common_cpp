// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
// clang-format on

template <typename T> std::pair<int, int> get_pair(int a, int b) {
  T a_T = UniversalScalarConversion<T, int>(a);
  T b_T = UniversalScalarConversion<T, int>(b);
  T res_T = ResInt(a_T, b_T);
  int res_int = UniversalScalarConversion<int, T>(res_T);
  T quo_T = QuoInt(a_T, b_T);
  int quo_int = UniversalScalarConversion<int, T>(quo_T);
  return {res_int, quo_int};
}

void check_consistency(
    int a, int b,
    std::vector<std::pair<std::pair<int, int>, std::string>> const &l_result,
    size_t &n_error) {
  size_t n_result = l_result.size();
  for (size_t i_res = 0; i_res < n_result; i_res++) {
    for (size_t j_res = i_res + 1; j_res < n_result; j_res++) {
      auto eP1 = l_result[i_res];
      auto eP2 = l_result[j_res];
      if (eP1.first.first != eP2.first.first ||
          eP1.first.second != eP2.first.second) {
        std::cerr << "Error for a=" << a << " b=" << b << "\n";
        std::cerr << "For class " << eP1.second
                  << " we found res=" << eP1.first.first
                  << " quot=" << eP1.first.second << "\n";
        n_error++;
      }
    }
  }
}

int main() {
  try {
    size_t n_error = 0;
    auto TestCons = [&](int a, int b) -> void {
      if (b != 0) {
        std::vector<std::pair<std::pair<int, int>, std::string>> l_result;
        l_result.push_back({get_pair<int64_t>(a, b), "int64_t"});
        l_result.push_back({get_pair<int32_t>(a, b), "int32_t"});
        l_result.push_back({get_pair<mpz_class>(a, b), "mpz_class"});
        l_result.push_back({get_pair<mpq_class>(a, b), "mpq_class"});
        l_result.push_back({get_pair<boost::multiprecision::cpp_int>(a, b),
                            "boost::multiprecision::cpp_int"});
        check_consistency(a, b, l_result, n_error);
      }
    };
    int nb = 100;
    int siz = 10000;
    for (int i = 0; i < nb; i++) {
      std::cerr << "i=" << i << "/" << nb << "\n";
      int a = random() % (2 * siz + 1) - siz;
      int b = random() % (2 * siz + 1) - siz;
      TestCons(a, b);
    }
    for (int a = -10; a < 10; a++)
      for (int b = -10; b < 10; b++)
        TestCons(a, b);
    std::cerr << "n_error=" << n_error << "\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
