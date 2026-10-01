// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYBOOSTCPPINT_H_
#define SRC_NUMBER_NUMBERTHEORYBOOSTCPPINT_H_
#define INCLUDE_NUMBER_THEORY_BOOST_CPP_INT
// clang-format off
#include "BasicNumberTypes.h"
#include "NumberTheoryBoostFormat.h"
#include "ExceptionsFunc.h"
#include "TemplateTraits.h"
#include "TypeConversion.h"
#include "boost_serialization.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/serialization/split_free.hpp>
#include <limits>
#include <string>
#include <utility>
// clang-format on

// The FMA form: cpp_int materializes a temporary for `a*b`, so the
// reused-scratch form is faster (measured).
// clang-format off
#define NT_BOOST_INT boost::multiprecision::cpp_int
#define NT_BOOST_RAT boost::multiprecision::cpp_rational
#define NT_BOOST_IS_INT is_boost_cpp_int
#define NT_BOOST_IS_RAT is_boost_cpp_rational
#define NT_BOOST_INT_FMA_PREFERED false
#define NT_BOOST_INT_TO_SMALL cpp_int_to_small_integer
#define NT_BOOST_CEIL_RAT Ceil_cpp_rational
#define NT_BOOST_FLOOR_RAT Floor_cpp_rational
#include "NumberTheoryBoost_impl.h"
// clang-format on

// No boost serialization here: unlike mpq_rational, cpp_rational and cpp_int
// have turned out not to need one.

// See T_frexp in TypeConversion.h. The 64 leading bits of |x| carry all that
// a double can hold; they are converted and the shift goes to the exponent.
inline double T_frexp(boost::multiprecision::cpp_int const &x, long &e) {
  if (x == 0) {
    e = 0;
    return 0.0;
  }
  boost::multiprecision::cpp_int x_abs = abs(x);
  long n_bit = static_cast<long>(boost::multiprecision::msb(x_abs)) + 1;
  long shift = n_bit > 64 ? n_bit - 64 : 0;
  boost::multiprecision::cpp_int top = x_abs >> shift;
  double m = static_cast<double>(top.convert_to<uint64_t>());
  if (x < 0) {
    m = -m;
  }
  return T_frexp_normalize(m, shift, e);
}

// See T_frexp in TypeConversion.h.
inline double T_frexp(boost::multiprecision::cpp_rational const &x, long &e) {
  long e_num, e_den;
  double m_num = T_frexp(boost::multiprecision::cpp_int(numerator(x)), e_num);
  double m_den =
      T_frexp(boost::multiprecision::cpp_int(denominator(x)), e_den);
  return T_frexp_quotient(m_num, e_num, m_den, e_den, e);
}

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYBOOSTCPPINT_H_
// clang-format on
