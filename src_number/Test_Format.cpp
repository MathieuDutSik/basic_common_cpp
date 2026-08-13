// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Exercises the std::formatter specializations of the number types.
//
// These replaced overloads of std::to_string that the project used to inject
// into namespace std. Adding a function to namespace std is undefined whatever
// its argument type ([namespace.std]/1), the permission of /2 covering
// specializations of standard class templates only. std::formatter is such a
// class template, so specializing it for a program-defined type is the
// sanctioned way of making the type printable.
//
// The point of this test is that the formatters be instantiated: an unused
// specialization compiles even when it is wrong, so every type is formatted
// here and the result compared against the stream output it must agree with.

// The boost headers come first on purpose. TypeConversionFinal.h guards the
// mixed gmp/boost conversions on both INCLUDE_NUMBER_THEORY_* macros being
// defined by the time it is first reached, and NumberTheory.h reaches it
// through NumberTheoryCommon.h. Same order as Test_Conversion.cpp.
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "NumberTheoryPadic.h"
#include "TypeConversionFinal.h"
#include "MAT_Matrix.h"
#include <format>
#include <iostream>
#include <sstream>
#include <string>
// clang-format on

static int n_fail = 0;

#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << " : " << #cond    \
                << "\n";                                                       \
      n_fail++;                                                                \
    }                                                                          \
  } while (0)

// std::format of a value has to agree with what operator<< writes.
template <typename T>
void check_agrees_with_stream(T const &val, std::string const &name) {
  std::ostringstream os;
  os << val;
  std::string via_stream = os.str();
  std::string via_format = std::format("{}", val);
  if (via_stream != via_format) {
    std::cerr << "FAIL " << name << ": stream gives '" << via_stream
              << "' but format gives '" << via_format << "'\n";
    n_fail++;
  }
}

template <typename T> void test_integer_type(std::string const &name) {
  check_agrees_with_stream(T(0), name);
  check_agrees_with_stream(T(1), name);
  check_agrees_with_stream(T(-17), name);
  T big(1);
  for (int i = 0; i < 40; i++) {
    big *= 10;
  }
  check_agrees_with_stream(big, name);
  // An unevaluated expression, not the number type itself. Both gmpxx and
  // boost multiprecision return expression templates from their arithmetic,
  // and std::format picks the formatter on the exact type without applying
  // the conversion to the number type that an overloaded function would get.
  // These two lines are the reason the formatters are declared over the
  // expression families rather than over the concrete types.
  check_agrees_with_stream(-big, name);
  check_agrees_with_stream(big + T(1), name);
  // Embedded in a larger format string, with several arguments.
  std::string s = std::format("[{} , {}]", T(3), T(-4));
  CHECK(s == "[3 , -4]");
}

template <typename Tnum, typename Tden>
void test_rational_type(std::string const &name) {
  check_agrees_with_stream(Tnum(0), name);
  check_agrees_with_stream(Tnum(Tden(22)) / Tnum(Tden(7)), name);
  check_agrees_with_stream(Tnum(Tden(-1)) / Tnum(Tden(3)), name);
  std::string s = std::format("q={}", Tnum(Tden(1)) / Tnum(Tden(2)));
  CHECK(s == "q=1/2");
}

int main() {
  test_integer_type<mpz_class>("mpz_class");
  test_integer_type<boost::multiprecision::cpp_int>("cpp_int");
  test_integer_type<boost::multiprecision::mpz_int>("mpz_int");

  test_rational_type<mpq_class, mpz_class>("mpq_class");
  test_rational_type<boost::multiprecision::cpp_rational,
                     boost::multiprecision::cpp_int>("cpp_rational");
  test_rational_type<boost::multiprecision::mpq_rational,
                     boost::multiprecision::mpz_int>("mpq_rational");

  // Rational<Tint> of the project.
  using RatInt = Rational<int64_t>;
  check_agrees_with_stream(RatInt(3), "Rational<int64_t>");
  CHECK(std::format("{}", RatInt(3)) == "3");

  // The formatters are reachable through a generic function too, which is the
  // case that matters for the templated code of the project.
  auto generic = []<typename T>(T const &x) { return std::format("<{}>", x); };
  CHECK(generic(mpz_class(5)) == "<5>");
  CHECK(generic(boost::multiprecision::cpp_int(6)) == "<6>");

  // Expressions, again through a generic function.
  mpz_class a(7), b(8);
  CHECK(std::format("{}", a * b) == "56");
  boost::multiprecision::cpp_int c(9), d(10);
  CHECK(std::format("{}", c * d) == "90");

  // The call sites in the library that format a value of a generic number
  // type. Each of them lives in a function template, so a broken formatter
  // there is not a compile error until something instantiates it -- which is
  // exactly how the Padic one got through. They are instantiated here.

  // MinMaxMatrix, in MAT_MatrixFund.h.
  MyMatrix<mpq_class> M(2, 2);
  M(0, 0) = mpq_class(1) / mpq_class(2);
  M(0, 1) = mpq_class(5);
  M(1, 0) = mpq_class(-3);
  M(1, 1) = mpq_class(0);
  CHECK(MinMaxMatrix(M) == "min/max=-3 / 5");

  // Padic_debug_print, in NumberTheoryPadic.h.
  Padic<mpz_class> pad =
      Padic_from_positive_integer(mpz_class(12), mpz_class(2));
  CHECK(std::format("{}", pad) ==
        std::format("eff_valuation={} precision={} coefficients={}",
                    pad.eff_valuation, pad.precision,
                    [&]() {
                      std::string s;
                      for (auto &v : pad.coefficients) {
                        s += " " + std::format("{}", v);
                      }
                      return s;
                    }()));
  std::ostringstream discard;
  Padic_debug_print(pad, discard);
  CHECK(discard.str() == std::format("{}", pad) + "\n");

  // TYPE_CONVERSION_IsInteger, in TypeConversionFinal.h: the message is built
  // only on the failing path.
  try {
    int64_t out;
    TYPE_CONVERSION_IsInteger(stc<mpq_class>{mpq_class(1) / mpq_class(3)}, out);
    CHECK(false && "TYPE_CONVERSION_IsInteger should have thrown");
  } catch (ConversionException const &e) {
    CHECK(e.val == "a1=1/3 is not an integer");
  }

  // TYPE_CONVERSION_Rational_T, in rational.h: likewise. Rational<SafeInt64>
  // and not Rational<int64_t>: the denominator is what gets formatted, and for
  // a builtin integer std::to_string would compile, so the instantiation would
  // not hold the formatter to anything.
  try {
    SafeInt64 out;
    TYPE_CONVERSION_Rational_T(
        stc<Rational<SafeInt64>>{Rational<SafeInt64>(1, 3)}, out);
    CHECK(false && "TYPE_CONVERSION_Rational_T should have thrown");
  } catch (ConversionException const &e) {
    CHECK(e.val == "The denominator should be 1. It is den = 3");
  }

  if (n_fail == 0) {
    std::cerr << "Test_Format: all checks passed.\n";
    return 0;
  }
  std::cerr << "Test_Format: " << n_fail << " check(s) failed.\n";
  return 1;
}
