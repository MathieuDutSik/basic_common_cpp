// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
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
        l_result.push_back({get_pair<SafeInt64>(a, b), "SafeInt64"});
        l_result.push_back({get_pair<mpz_class>(a, b), "mpz_class"});
        l_result.push_back({get_pair<mpq_class>(a, b), "mpq_class"});
        l_result.push_back({get_pair<Rational<SafeInt64>>(a, b),
                            "Rational<SafeInt64>"});
        l_result.push_back({get_pair<boost::multiprecision::cpp_int>(a, b),
                            "boost::multiprecision::cpp_int"});
        l_result.push_back({get_pair<boost::multiprecision::cpp_rational>(a, b),
                            "boost::multiprecision::cpp_rational"});
        l_result.push_back({get_pair<boost::multiprecision::mpz_int>(a, b),
                            "boost::multiprecision::mpz_int"});
        l_result.push_back({get_pair<boost::multiprecision::mpq_rational>(a, b),
                            "boost::multiprecision::mpq_rational"});
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
    // Unsigned ResInt. QUO_INT is not defined for unsigned types, so only
    // ResInt is exercised, over non-negative arguments, across the unsigned
    // widths and against the reference a % b. This covers in particular the
    // uint64_t kernel, including arguments beyond 2^32 where the uint32_t
    // path cannot reach.
    auto TestUnsignedResInt = [&](uint64_t a, uint64_t b) -> void {
      if (b == 0)
        return;
      uint64_t ref = a % b;
      uint64_t r64 = ResInt<uint64_t>(a, b);
      if (r64 != ref) {
        std::cerr << "uint64_t ResInt error a=" << a << " b=" << b
                  << " got=" << r64 << " ref=" << ref << "\n";
        n_error++;
      }
      if (a <= std::numeric_limits<uint32_t>::max() &&
          b <= std::numeric_limits<uint32_t>::max()) {
        uint32_t r32 = ResInt<uint32_t>(static_cast<uint32_t>(a),
                                        static_cast<uint32_t>(b));
        if (static_cast<uint64_t>(r32) != ref) {
          std::cerr << "uint32_t ResInt error a=" << a << " b=" << b << "\n";
          n_error++;
        }
      }
    };
    for (uint64_t a = 0; a < 50; a++)
      for (uint64_t b = 1; b < 50; b++)
        TestUnsignedResInt(a, b);
    std::vector<uint64_t> big_vals = {45580ULL,
                                      91160ULL,
                                      4294967311ULL,
                                      1000000000000ULL,
                                      18446744073709551557ULL};
    for (uint64_t a : big_vals)
      for (uint64_t b : big_vals)
        TestUnsignedResInt(a, b);
    std::cerr << "n_error=" << n_error << "\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
