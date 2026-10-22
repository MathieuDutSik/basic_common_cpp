// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_ARITHMETIC_DISPATCH_TYPES_H_
#define SRC_NUMBER_ARITHMETIC_DISPATCH_TYPES_H_

// clang-format off
#include "Arithmetic_dispatch.h"
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#ifdef QUADRATIC_FIELDS
#include "NumberTheoryQuadField.h"
#endif
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "Basic_file.h"
#include "Basic_string.h"
#include <optional>
#include <string>
// clang-format on

/*
  The registration in an ArithmeticDispatch of the common families of
  arithmetics. The argument f is a generic lambda run for the chosen type:
    auto f = [&]<typename T>() -> void { process<T>(...); };
  Each optional type is registered under the macro that compiles it
  (ENABLE_BOOST_TYPES, ENABLE_FLINT_SUPPORT, QUADRATIC_FIELDS), so that the
  printed list is the compiled one.
*/

// The rational fields: mpq_class, the boost ones, fmpq_class. The type run
// for the name "mpq_class" is Tmpq, which a program can set to another
// implementation of the GMP rationals (e.g. the boost bindings for Oscar).
template <typename Tmpq = mpq_class, typename F>
void add_rational_arithmetics(ArithmeticDispatch &dispatch, F f) {
  dispatch.add("mpq_class", "rational arithmetic over GMP mpq_class",
               [f]() -> void { f.template operator()<Tmpq>(); });
#ifdef ENABLE_BOOST_TYPES
  dispatch.add("mpq_rational",
               "rational arithmetic over boost::multiprecision::mpq_rational",
               [f]() -> void {
                 f.template operator()<boost::multiprecision::mpq_rational>();
               });
  dispatch.add("cpp_rational",
               "rational arithmetic over boost::multiprecision::cpp_rational",
               [f]() -> void {
                 f.template operator()<boost::multiprecision::cpp_rational>();
               });
#endif
#ifdef ENABLE_FLINT_SUPPORT
  dispatch.add("fmpq_class", "rational arithmetic over FLINT fmpq_class",
               [f]() -> void { f.template operator()<fmpq_class>(); });
#endif
}

// The same rational fields, under the names of their libraries, as used by
// a number of programs: gmp, gmp_boost, multi_boost, flint.
template <typename F>
void add_rational_arithmetics_by_library(ArithmeticDispatch &dispatch, F f) {
  dispatch.add("gmp", "rational arithmetic over GMP mpq_class",
               [f]() -> void { f.template operator()<mpq_class>(); });
#ifdef ENABLE_BOOST_TYPES
  dispatch.add("gmp_boost",
               "rational arithmetic over boost::multiprecision::mpq_rational",
               [f]() -> void {
                 f.template operator()<boost::multiprecision::mpq_rational>();
               });
  dispatch.add("multi_boost",
               "rational arithmetic over boost::multiprecision::cpp_rational",
               [f]() -> void {
                 f.template operator()<boost::multiprecision::cpp_rational>();
               });
#endif
#ifdef ENABLE_FLINT_SUPPORT
  dispatch.add("flint", "rational arithmetic over FLINT fmpq_class",
               [f]() -> void { f.template operator()<fmpq_class>(); });
#endif
}

// The rational fields with their integer rings, under the names of their
// libraries; f is run as f.template operator()<T, Tint>(). A program whose
// code does not support the boost types passes with_boost = false, so that a
// build with ENABLE_BOOST_TYPES does not instantiate it on them.
template <bool with_boost = true, typename F>
void add_rational_integer_arithmetics_by_library(ArithmeticDispatch &dispatch,
                                                 F f) {
  dispatch.add("gmp", "GMP mpq_class and mpz_class", [f]() -> void {
    f.template operator()<mpq_class, mpz_class>();
  });
#ifdef ENABLE_BOOST_TYPES
  if constexpr (with_boost) {
    dispatch.add("gmp_boost",
                 "boost::multiprecision::mpq_rational and mpz_int",
                 [f]() -> void {
                   f.template operator()<boost::multiprecision::mpq_rational,
                                         boost::multiprecision::mpz_int>();
                 });
    dispatch.add("multi_boost",
                 "boost::multiprecision::cpp_rational and cpp_int",
                 [f]() -> void {
                   f.template operator()<boost::multiprecision::cpp_rational,
                                         boost::multiprecision::cpp_int>();
                 });
  }
#endif
#ifdef ENABLE_FLINT_SUPPORT
  dispatch.add("flint", "FLINT fmpq_class and fmpz_class", [f]() -> void {
    f.template operator()<fmpq_class, fmpz_class>();
  });
#endif
}

// The namelist documentation of the arithmetic key of the programs using
// add_rational_integer_arithmetics_by_library: flint by default when
// compiled in, as it is usually faster, gmp otherwise.
template <bool with_boost = true>
std::string doc_rational_integer_arithmetics_by_library() {
  ArithmeticDispatch dispatch;
  add_rational_integer_arithmetics_by_library<with_boost>(
      dispatch, []<typename T, typename Tint>() -> void {});
  return dispatch.doc({"flint", "gmp"});
}

// The integer rings: mpz_class, the boost ones, fmpz_class.
template <typename F>
void add_integer_arithmetics(ArithmeticDispatch &dispatch, F f) {
  dispatch.add("mpz_class", "integer arithmetic over GMP mpz_class",
               [f]() -> void { f.template operator()<mpz_class>(); });
#ifdef ENABLE_BOOST_TYPES
  dispatch.add("mpz_int",
               "integer arithmetic over boost::multiprecision::mpz_int",
               [f]() -> void {
                 f.template operator()<boost::multiprecision::mpz_int>();
               });
  dispatch.add("cpp_int",
               "integer arithmetic over boost::multiprecision::cpp_int",
               [f]() -> void {
                 f.template operator()<boost::multiprecision::cpp_int>();
               });
#endif
#ifdef ENABLE_FLINT_SUPPORT
  dispatch.add("fmpz_class", "integer arithmetic over FLINT fmpz_class",
               [f]() -> void { f.template operator()<fmpz_class>(); });
#endif
}

// The algebraic fields: the quadratic fields Q(sqrt(2)) and Q(sqrt(5)), and
// the real algebraic field described in a file.
template <typename F>
void add_algebraic_field_arithmetics(ArithmeticDispatch &dispatch, F f) {
#ifdef QUADRATIC_FIELDS
  dispatch.add("Qsqrt2", "arithmetic over the field Q(sqrt(2))",
               [f]() -> void {
                 f.template operator()<QuadField<Trat_quad_field, 2>>();
               });
  dispatch.add("Qsqrt5", "arithmetic over the field Q(sqrt(5))",
               [f]() -> void {
                 f.template operator()<QuadField<Trat_quad_field, 5>>();
               });
#endif
  dispatch.add_prefix(
      "RealAlgebraic=", "FileDesc",
      "arithmetic over the real algebraic field described in FileDesc",
      [f](std::string const &FileAlgebraicField) -> void {
        if (!FILE_IsExistingFile(FileAlgebraicField)) {
          std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                    << " is missing\n";
          throw TerminalException{1};
        }
        HelperClassRealField<Trat_real_field> hcrf(FileAlgebraicField);
        int const idx_real_algebraic_field = 1;
        insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
        f.template operator()<RealField<idx_real_algebraic_field>>();
      });
}

// clang-format off
#endif  // SRC_NUMBER_ARITHMETIC_DISPATCH_TYPES_H_
// clang-format on
