// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MATRIXTYPES_H_
#define SRC_MATRIX_MATRIXTYPES_H_

// clang-format off
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/Sparse>
// clang-format on

#if defined INCLUDE_NUMBER_THEORY_BOOST_GMP_INT ||                             \
    defined INCLUDE_NUMBER_THEORY_BOOST_CPP_INT
#include <boost/multiprecision/eigen.hpp>
#endif

template <typename T> using MyVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename T>
using MyMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
using MySparseMatrix = Eigen::SparseMatrix<T, Eigen::ColMajor>;

// clang-format off
#endif  // SRC_MATRIX_MATRIXTYPES_H_
// clang-format on
