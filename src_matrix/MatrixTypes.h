// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MATRIXTYPES_H_
#define SRC_MATRIX_MATRIXTYPES_H_

// clang-format off
#include <Eigen/Dense>
#include <Eigen/LU>
// clang-format on

#if defined INCLUDE_NUMBER_THEORY_BOOST_GMP_INT ||                             \
    defined INCLUDE_NUMBER_THEORY_BOOST_CPP_INT
#include <boost/multiprecision/eigen.hpp>
#endif

#if defined INCLUDE_NUMBER_THEORY_BOOST_CPP_INT
// Must come after boost/multiprecision/eigen.hpp: overrides Literal=double
// for cpp_int / cpp_rational so Eigen >= 3.5 does not emit
// (cpp_int == double) inside is_exactly_zero / equal_strict.
#include "EigenBoostNumTraits.h"
#endif

template <typename T> using MyVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename T>
using MyMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

// clang-format off
#endif  // SRC_MATRIX_MATRIXTYPES_H_
// clang-format on
