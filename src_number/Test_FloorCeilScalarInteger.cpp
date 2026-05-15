// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "TypeConversion.h"
#include "TypeConversionFinal.h"
// clang-format on

// Verify, for a rational type Trat and an integer type Tint, that
// UniversalFloorScalarInteger and UniversalCeilScalarInteger satisfy the
// expected mathematical properties on a fixed sample of inputs.
template <typename Trat, typename Tint>
void check_pair(std::string const &name) {
  // Sample: numerator, denominator, expected floor, expected ceil.
  struct Case {
    int num, den, floor_expected, ceil_expected;
  };
  std::vector<Case> cases = {
      {7, 3, 2, 3},     {-7, 3, -3, -2}, {5, 2, 2, 3},
      {-5, 2, -3, -2},  {4, 2, 2, 2},    {-4, 2, -2, -2},
      {0, 1, 0, 0},     {1, 1, 1, 1},    {-1, 1, -1, -1},
      {100, 7, 14, 15}, {-100, 7, -15, -14},
  };
  for (auto const &c : cases) {
    Trat num(c.num);
    Trat den(c.den);
    Trat x = num / den;
    Tint got_floor = UniversalFloorScalarInteger<Tint, Trat>(x);
    Tint got_ceil = UniversalCeilScalarInteger<Tint, Trat>(x);
    Tint exp_floor(c.floor_expected);
    Tint exp_ceil(c.ceil_expected);
    if (got_floor != exp_floor) {
      std::cerr << name << ": floor(" << x << ") gave " << got_floor
                << " expected " << exp_floor << "\n";
      throw TerminalException{1};
    }
    if (got_ceil != exp_ceil) {
      std::cerr << name << ": ceil(" << x << ") gave " << got_ceil
                << " expected " << exp_ceil << "\n";
      throw TerminalException{1};
    }
    // Sanity: floor <= x <= ceil and ceil - floor in {0, 1}.
    if (Trat(got_floor) > x || Trat(got_ceil) < x) {
      std::cerr << name << ": floor/ceil bracketing failed for x=" << x
                << "\n";
      throw TerminalException{1};
    }
    Tint diff = got_ceil - got_floor;
    if (diff != Tint(0) && diff != Tint(1)) {
      std::cerr << name << ": ceil - floor not in {0,1} for x=" << x << "\n";
      throw TerminalException{1};
    }
  }
  std::cerr << "OK: " << name << "\n";
}

int main() {
  try {
    // mpq_class: floor/ceil to mpz_class, int, long
    check_pair<mpq_class, mpz_class>("mpq_class -> mpz_class");
    check_pair<mpq_class, int>("mpq_class -> int");
    check_pair<mpq_class, long>("mpq_class -> long");

    // Rational<SafeInt64>: floor/ceil to SafeInt64, int, long
    check_pair<Rational<SafeInt64>, SafeInt64>(
        "Rational<SafeInt64> -> SafeInt64");
    check_pair<Rational<SafeInt64>, int>("Rational<SafeInt64> -> int");
    check_pair<Rational<SafeInt64>, long>("Rational<SafeInt64> -> long");

    // boost::cpp_rational: floor/ceil to cpp_int, int, long
    check_pair<boost::multiprecision::cpp_rational,
               boost::multiprecision::cpp_int>("cpp_rational -> cpp_int");
    check_pair<boost::multiprecision::cpp_rational, int>(
        "cpp_rational -> int");
    check_pair<boost::multiprecision::cpp_rational, long>(
        "cpp_rational -> long");

    // boost::mpq_rational: floor/ceil to mpz_int, int, long
    check_pair<boost::multiprecision::mpq_rational,
               boost::multiprecision::mpz_int>("mpq_rational -> mpz_int");
    check_pair<boost::multiprecision::mpq_rational, int>(
        "mpq_rational -> int");
    check_pair<boost::multiprecision::mpq_rational, long>(
        "mpq_rational -> long");

    std::cerr << "Normal termination of Test_FloorCeilScalarInteger\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Test_FloorCeilScalarInteger\n";
    exit(e.eVal);
  }
}
