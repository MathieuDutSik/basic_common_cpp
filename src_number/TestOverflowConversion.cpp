// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Checks that UniversalScalarConversionCheck reports the overflow instead of
// returning a truncated value, for every multiprecision source type and every
// bounded integer target.
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheorySafeInt.h"
#include <limits>
#include <string>
// clang-format on

static int n_error = 0;

template <typename Tout, typename Tin>
void check_one(std::string const &name_in, std::string const &name_out) {
  std::string desc = name_in + " -> " + name_out;
  // The values that fit have to be converted, and to round-trip.
  std::vector<int64_t> l_inside{0, 1, -1};
  int64_t out_min = static_cast<int64_t>(std::numeric_limits<Tout>::min());
  int64_t out_max = static_cast<int64_t>(std::numeric_limits<Tout>::max());
  if constexpr (std::is_unsigned_v<Tout>) {
    // max may not fit in int64_t for uint64_t, handled by the outside part.
    if (sizeof(Tout) < 8) {
      l_inside.push_back(out_max);
    }
    l_inside.erase(l_inside.begin() + 2);
  } else {
    l_inside.push_back(out_min);
    l_inside.push_back(out_max);
  }
  for (auto &val : l_inside) {
    Tin val_in = UniversalScalarConversion<Tin, int64_t>(val);
    std::optional<Tout> opt = UniversalScalarConversionCheck<Tout, Tin>(val_in);
    if (!opt) {
      std::cerr << "ERROR " << desc << ": value " << val
                << " fits but was rejected\n";
      n_error++;
      continue;
    }
    int64_t back = static_cast<int64_t>(*opt);
    if (back != val) {
      std::cerr << "ERROR " << desc << ": value " << val << " came back as "
                << back << "\n";
      n_error++;
    }
  }
  // The values that do not fit have to be rejected.
  std::vector<Tin> l_outside;
  Tin one = UniversalScalarConversion<Tin, int64_t>(1);
  Tin big = UniversalScalarConversion<Tin, int64_t>(1);
  for (int i = 0; i < 200; i++) {
    big = big * 2;
  }
  l_outside.push_back(big);
  l_outside.push_back(-big);
  if constexpr (std::is_unsigned_v<Tout>) {
    l_outside.push_back(-one);
    if (sizeof(Tout) < 8) {
      l_outside.push_back(
          UniversalScalarConversion<Tin, int64_t>(out_max) + one);
    }
  } else {
    l_outside.push_back(UniversalScalarConversion<Tin, int64_t>(out_max) + one);
    l_outside.push_back(UniversalScalarConversion<Tin, int64_t>(out_min) - one);
  }
  for (auto &val_in : l_outside) {
    std::optional<Tout> opt = UniversalScalarConversionCheck<Tout, Tin>(val_in);
    if (opt) {
      std::cerr << "ERROR " << desc << ": out of range value was accepted as "
                << static_cast<int64_t>(*opt) << "\n";
      n_error++;
    }
  }
}

template <typename Tin> void check_all_targets(std::string const &name_in) {
  check_one<int8_t, Tin>(name_in, "int8_t");
  check_one<uint8_t, Tin>(name_in, "uint8_t");
  check_one<int16_t, Tin>(name_in, "int16_t");
  check_one<uint16_t, Tin>(name_in, "uint16_t");
  check_one<int32_t, Tin>(name_in, "int32_t");
  check_one<uint32_t, Tin>(name_in, "uint32_t");
  check_one<int64_t, Tin>(name_in, "int64_t");
  check_one<uint64_t, Tin>(name_in, "uint64_t");
}

int main() {
  try {
    check_all_targets<mpz_class>("mpz_class");
    check_all_targets<mpq_class>("mpq_class");
    check_all_targets<boost::multiprecision::mpz_int>("boost mpz_int");
    check_all_targets<boost::multiprecision::mpq_rational>("boost mpq");
    check_all_targets<boost::multiprecision::cpp_int>("boost cpp_int");
    check_all_targets<boost::multiprecision::cpp_rational>("boost cpp_rat");
    // A non integral rational has to be rejected as well.
    mpq_class val_frac(1, 2);
    if (UniversalScalarConversionCheck<int64_t, mpq_class>(val_frac)) {
      std::cerr << "ERROR mpq_class -> int64_t: 1/2 was accepted\n";
      n_error++;
    }
    if (n_error == 0) {
      std::cerr << "ALL CONVERSION TESTS PASSED\n";
      return 0;
    }
    std::cerr << "n_error=" << n_error << "\n";
    return 1;
  } catch (TerminalException const &e) {
    std::cerr << "Error in TestOverflowConversion\n";
    exit(e.eVal);
  }
}
