// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MATRIXTYPESSPARSE_H_
#define SRC_MATRIX_MATRIXTYPESSPARSE_H_

#include <Eigen/Sparse>

#if defined INCLUDE_NUMBER_THEORY_BOOST_GMP_INT ||                             \
    defined INCLUDE_NUMBER_THEORY_BOOST_CPP_INT
#include <boost/multiprecision/eigen.hpp>
#endif

template <typename T>
using MySparseMatrix = Eigen::SparseMatrix<T, Eigen::ColMajor>;

// clang-format off
#endif  // SRC_MATRIX_MATRIXTYPESSPARSE_H_
// clang-format on
