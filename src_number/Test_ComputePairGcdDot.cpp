// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheorySafeInt.h"
#include "TypeConversion.h"
#include <random>
// clang-format on

template <typename T> void Test_ComputePairGcdDot() {
  std::random_device rd;
  std::mt19937 gen(rd());

  // Test with various ranges
  std::vector<std::pair<int, int>> ranges = {
    {-100, 100},
    {-1000, 1000},
    {-10000, 10000}
  };

  size_t n_error = 0;
  size_t n_tests = 0;

  for (auto const& range : ranges) {
    std::uniform_int_distribution<int> dis(range.first, range.second);

    for (int i = 0; i < 100; i++) {
      int m_int = dis(gen);
      int n_int = dis(gen);

      // Skip the case where both are zero (handled separately)
      if (m_int == 0 && n_int == 0) continue;

      T m = UniversalScalarConversion<T, int>(m_int);
      T n = UniversalScalarConversion<T, int>(n_int);

      PairGCD_dot<T> result = ComputePairGcdDot(m, n);
      T eCoeff1 = result.a;
      T eCoeff2 = result.b;
      T gcd = result.gcd;

      n_tests++;

      // Check the fundamental property: eCoeff1 * m + eCoeff2 * n = gcd
      T check_gcd = eCoeff1 * m + eCoeff2 * n;
      std::cerr << "m=" << m << " n=" << n << " eCoeff12=" << eCoeff1 << " / " << eCoeff2 << "\n";
      if (check_gcd != gcd) {
        std::cerr << "ERROR: Extended Euclidean property failed\n";
        std::cerr << "m=" << m << ", n=" << n << "\n";
        std::cerr << "eCoeff1=" << eCoeff1 << ", eCoeff2=" << eCoeff2 << "\n";
        std::cerr << "gcd=" << gcd << ", computed=" << check_gcd << "\n";
        n_error++;
        continue;
      }

      // Check that coefficients are "small enough"
      // For Extended Euclidean Algorithm, coefficients should satisfy |eCoeff1| <= |n/gcd|
      // and |eCoeff2| <= |m/gcd| when both inputs are non-zero
      if (m != 0 && n != 0 && gcd != 0) {
        T abs_m = (m < 0) ? -m : m;
        T abs_n = (n < 0) ? -n : n;
        T abs_eCoeff1 = (eCoeff1 < 0) ? -eCoeff1 : eCoeff1;
        T abs_eCoeff2 = (eCoeff2 < 0) ? -eCoeff2 : eCoeff2;
        T abs_gcd = (gcd < 0) ? -gcd : gcd;

        T bound1 = abs_n / abs_gcd;
        T bound2 = abs_m / abs_gcd;

        if (abs_eCoeff1 > bound1) {
          std::cerr << "WARNING: |eCoeff1| = " << abs_eCoeff1 << " > " << bound1 << " = |n|/|gcd|\n";
          std::cerr << "m=" << m << ", n=" << n << ", gcd=" << gcd << "\n";
          n_error++;
        }

        if (abs_eCoeff2 > bound2) {
          std::cerr << "WARNING: |eCoeff2| = " << abs_eCoeff2 << " > " << bound2 << " = |m|/|gcd|\n";
          std::cerr << "m=" << m << ", n=" << n << ", gcd=" << gcd << "\n";
          n_error++;
        }
      }
    }
  }

  // Test edge cases
  T zero = UniversalScalarConversion<T, int>(0);
  T one = UniversalScalarConversion<T, int>(1);
  T minus_one = UniversalScalarConversion<T, int>(-1);

  std::vector<std::pair<T, T>> edge_cases = {
    {one, zero}, {zero, one}, {minus_one, zero}, {zero, minus_one},
    {one, one}, {one, minus_one}, {minus_one, one}, {minus_one, minus_one}
  };

  for (auto const& edge_case : edge_cases) {
    T m = edge_case.first;
    T n = edge_case.second;

    PairGCD_dot<T> result = ComputePairGcdDot(m, n);
    T eCoeff1 = result.a;
    T eCoeff2 = result.b;
    T gcd = result.gcd;

    n_tests++;

    // Check the fundamental property
    T check_gcd = eCoeff1 * m + eCoeff2 * n;
    if (check_gcd != gcd) {
      std::cerr << "ERROR: Extended Euclidean property failed for edge case\n";
      std::cerr << "m=" << m << ", n=" << n << "\n";
      std::cerr << "eCoeff1=" << eCoeff1 << ", eCoeff2=" << eCoeff2 << "\n";
      std::cerr << "gcd=" << gcd << ", computed=" << check_gcd << "\n";
      n_error++;
    }
  }

  std::cerr << "Test completed for numerical type with " << n_tests << " tests\n";
  std::cerr << "Number of errors: " << n_error << "\n";

  if (n_error > 0) {
    throw TerminalException{1};
  }
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 2) {
      std::cerr << "Test_ComputePairGcdDot [numeric_type]\n";
      std::cerr << "\n";
      std::cerr << "numeric_type options:\n";
      std::cerr << "  int64_t\n";
      std::cerr << "  mpz_class\n";
      std::cerr << "  SafeInt64\n";
      std::cerr << "  boost_cpp_int\n";
      throw TerminalException{1};
    }

    std::string numeric_type = argv[1];

    std::cerr << "Testing ComputePairGcdDot with type: " << numeric_type << "\n";

    if (numeric_type == "int64_t") {
      Test_ComputePairGcdDot<int64_t>();
    } else if (numeric_type == "mpz_class") {
      Test_ComputePairGcdDot<mpz_class>();
    } else if (numeric_type == "SafeInt64") {
      Test_ComputePairGcdDot<SafeInt64>();
    } else if (numeric_type == "boost_cpp_int") {
      Test_ComputePairGcdDot<boost::multiprecision::cpp_int>();
    } else {
      std::cerr << "Unknown numeric type: " << numeric_type << "\n";
      throw TerminalException{1};
    }

    std::cerr << "All tests passed successfully!\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
