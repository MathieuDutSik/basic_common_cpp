// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
// clang-format on

struct OpResults {
  int sum;   // a + b
  int prod;  // a * b
  int inc;   // a += 1
  int mul2;  // a *= 2
  int sub;   // a -= c
};

template <typename T> OpResults get_ops(int a, int b, int c) {
  T a_T = UniversalScalarConversion<T, int>(a);
  T b_T = UniversalScalarConversion<T, int>(b);
  T c_T = UniversalScalarConversion<T, int>(c);
  T sum_T = a_T + b_T;
  T prod_T = a_T * b_T;
  T inc_T = a_T;
  inc_T += 1;
  T mul2_T = a_T;
  mul2_T *= 2;
  T sub_T = a_T;
  sub_T -= c_T;
  OpResults res;
  res.sum = UniversalScalarConversion<int, T>(sum_T);
  res.prod = UniversalScalarConversion<int, T>(prod_T);
  res.inc = UniversalScalarConversion<int, T>(inc_T);
  res.mul2 = UniversalScalarConversion<int, T>(mul2_T);
  res.sub = UniversalScalarConversion<int, T>(sub_T);
  return res;
}

bool eq_ops(OpResults const &x, OpResults const &y) {
  return x.sum == y.sum && x.prod == y.prod && x.inc == y.inc &&
         x.mul2 == y.mul2 && x.sub == y.sub;
}

void print_ops(std::ostream &os, OpResults const &r) {
  os << "sum=" << r.sum << " prod=" << r.prod << " inc=" << r.inc
     << " mul2=" << r.mul2 << " sub=" << r.sub;
}

void check_consistency(
    int a, int b, int c,
    std::vector<std::pair<OpResults, std::string>> const &l_result,
    size_t &n_error) {
  size_t n_result = l_result.size();
  for (size_t i_res = 0; i_res < n_result; i_res++) {
    for (size_t j_res = i_res + 1; j_res < n_result; j_res++) {
      auto const &eP1 = l_result[i_res];
      auto const &eP2 = l_result[j_res];
      if (!eq_ops(eP1.first, eP2.first)) {
        std::cerr << "Error for a=" << a << " b=" << b << " c=" << c << "\n";
        std::cerr << "  " << eP1.second << ": ";
        print_ops(std::cerr, eP1.first);
        std::cerr << "\n";
        std::cerr << "  " << eP2.second << ": ";
        print_ops(std::cerr, eP2.first);
        std::cerr << "\n";
        n_error++;
      }
    }
  }
}

int main() {
  try {
    size_t n_error = 0;
    auto TestCons = [&](int a, int b, int c) -> void {
      std::vector<std::pair<OpResults, std::string>> l_result;
      l_result.push_back({get_ops<int64_t>(a, b, c), "int64_t"});
      l_result.push_back({get_ops<int32_t>(a, b, c), "int32_t"});
      l_result.push_back({get_ops<SafeInt64>(a, b, c), "SafeInt64"});
      l_result.push_back({get_ops<mpz_class>(a, b, c), "mpz_class"});
      l_result.push_back({get_ops<mpq_class>(a, b, c), "mpq_class"});
      l_result.push_back(
          {get_ops<Rational<SafeInt64>>(a, b, c), "Rational<SafeInt64>"});
      l_result.push_back({get_ops<boost::multiprecision::cpp_int>(a, b, c),
                          "boost::multiprecision::cpp_int"});
      l_result.push_back({get_ops<boost::multiprecision::cpp_rational>(a, b, c),
                          "boost::multiprecision::cpp_rational"});
      l_result.push_back({get_ops<boost::multiprecision::mpz_int>(a, b, c),
                          "boost::multiprecision::mpz_int"});
      l_result.push_back(
          {get_ops<boost::multiprecision::mpq_rational>(a, b, c),
           "boost::multiprecision::mpq_rational"});
      check_consistency(a, b, c, l_result, n_error);
    };
    // Bound chosen so that a*b fits comfortably in int32 (worst case 10000).
    int nb = 100;
    int siz = 100;
    for (int i = 0; i < nb; i++) {
      int a = random() % (2 * siz + 1) - siz;
      int b = random() % (2 * siz + 1) - siz;
      int c = random() % (2 * siz + 1) - siz;
      TestCons(a, b, c);
    }
    for (int a = -10; a <= 10; a++)
      for (int b = -10; b <= 10; b++)
        for (int c = -10; c <= 10; c++)
          TestCons(a, b, c);
    std::cerr << "n_error=" << n_error << "\n";
    if (n_error != 0) {
      throw TerminalException{1};
    }
    std::cerr << "Normal termination of Test_Operations\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Test_Operations\n";
    exit(e.eVal);
  }
}
