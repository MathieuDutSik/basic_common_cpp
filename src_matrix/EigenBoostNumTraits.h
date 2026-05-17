// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_EIGENBOOSTNUMTRAITS_H_
#define SRC_MATRIX_EIGENBOOSTNUMTRAITS_H_

// Override Eigen::NumTraits<...>::Literal for boost::multiprecision integer
// and rational backends.
//
// Background: boost/multiprecision/eigen.hpp sets Literal = double for every
// number<Backend, ET>. Starting with Eigen 3.5 (= upstream "Eigen 5.0"),
// is_exactly_zero(x) is implemented as equal_strict(x, NumTraits<X>::Literal{0})
// and equal_strict ultimately writes `return x == y;`. For exact integer /
// rational types this expands to `cpp_int == double{0}`, which boost
// cpp_int / cpp_rational rightly do not provide an operator== for.
//
// We retarget Literal to `int`, which:
//   - matches what Eigen does for built-in integer scalars
//     (NumTraits<int>::Literal == int, NumTraits<long>::Literal == long);
//   - costs nothing (no number<Backend> is constructed for `Literal{0}`);
//   - relies on boost's existing operator==(number<Backend, ET>, int),
//     never on cpp_int <-> double comparisons.
//
// Crucially, this does NOT introduce any new operator: user-side code such
// as `cpp_int x; bool b = (x == 0.0);` continues to fail to compile, which
// is exactly what one wants from exact-arithmetic types.
//
// Coverage is opt-in via the same INCLUDE_NUMBER_THEORY_BOOST_* macros that
// gate boost/multiprecision/eigen.hpp itself, so this header never drags in
// GMP unless the caller already asked for it.

#include <Eigen/Core>
#include <boost/multiprecision/eigen.hpp>

#if defined INCLUDE_NUMBER_THEORY_BOOST_CPP_INT
#include <boost/multiprecision/cpp_int.hpp>

namespace Eigen {

template <boost::multiprecision::expression_template_option ET>
struct NumTraits<boost::multiprecision::number<
    boost::multiprecision::backends::cpp_int_backend<>, ET>>
    : NumTraitsImp<boost::multiprecision::number<
                       boost::multiprecision::backends::cpp_int_backend<>, ET>,
                   typename boost::multiprecision::number<
                       boost::multiprecision::backends::cpp_int_backend<>,
                       ET>::value_type> {
  using Literal = int;
};

template <boost::multiprecision::expression_template_option ET>
struct NumTraits<boost::multiprecision::number<
    boost::multiprecision::backends::rational_adaptor<
        boost::multiprecision::backends::cpp_int_backend<>>,
    ET>>
    : NumTraitsImp<boost::multiprecision::number<
                       boost::multiprecision::backends::rational_adaptor<
                           boost::multiprecision::backends::cpp_int_backend<>>,
                       ET>,
                   typename boost::multiprecision::number<
                       boost::multiprecision::backends::rational_adaptor<
                           boost::multiprecision::backends::cpp_int_backend<>>,
                       ET>::value_type> {
  using Literal = int;
};

}  // namespace Eigen
#endif  // INCLUDE_NUMBER_THEORY_BOOST_CPP_INT

#if defined INCLUDE_NUMBER_THEORY_BOOST_GMP_INT
#include <boost/multiprecision/gmp.hpp>

namespace Eigen {

template <boost::multiprecision::expression_template_option ET>
struct NumTraits<
    boost::multiprecision::number<boost::multiprecision::backends::gmp_int, ET>>
    : NumTraitsImp<boost::multiprecision::number<
                       boost::multiprecision::backends::gmp_int, ET>,
                   typename boost::multiprecision::number<
                       boost::multiprecision::backends::gmp_int,
                       ET>::value_type> {
  using Literal = int;
};

template <boost::multiprecision::expression_template_option ET>
struct NumTraits<boost::multiprecision::number<
    boost::multiprecision::backends::gmp_rational, ET>>
    : NumTraitsImp<boost::multiprecision::number<
                       boost::multiprecision::backends::gmp_rational, ET>,
                   typename boost::multiprecision::number<
                       boost::multiprecision::backends::gmp_rational,
                       ET>::value_type> {
  using Literal = int;
};

}  // namespace Eigen
#endif  // INCLUDE_NUMBER_THEORY_BOOST_GMP_INT

// clang-format off
#endif  // SRC_MATRIX_EIGENBOOSTNUMTRAITS_H_
// clang-format on
