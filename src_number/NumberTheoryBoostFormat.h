// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYBOOSTFORMAT_H_
#define SRC_NUMBER_NUMBERTHEORYBOOSTFORMAT_H_
// clang-format off
#include "BasicNumberTypes.h"
#include <boost/multiprecision/number.hpp>
// clang-format on

/*
  std::format support for the boost multiprecision types.

  Two specializations rather than one per concrete type. The first covers every
  number<Backend, ET>, so cpp_int, cpp_rational, mpz_int and mpq_rational at
  once. The second covers the expression templates: the arithmetic of boost
  multiprecision does not return the number type but an
  expression<tag, ...> that converts to it, exactly as gmpxx does, and the
  formatter of std::format is picked on the exact type with no conversion
  allowed. Without the second one, std::format("{}", a + b) would not compile
  even though std::format("{}", a) does.

  They live in a header of their own because both NumberTheoryBoostCppInt.h and
  NumberTheoryBoostGmpInt.h need them and the two can be included together; a
  copy in each would be a redefinition.
 */
template <typename Backend, boost::multiprecision::expression_template_option ET>
struct std::formatter<boost::multiprecision::number<Backend, ET>>
    : ostream_formatter<boost::multiprecision::number<Backend, ET>> {};

template <typename tag, typename A1, typename A2, typename A3, typename A4>
struct std::formatter<boost::multiprecision::detail::expression<tag, A1, A2, A3, A4>>
    : ostream_formatter<
          boost::multiprecision::detail::expression<tag, A1, A2, A3, A4>> {};

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYBOOSTFORMAT_H_
// clang-format on
