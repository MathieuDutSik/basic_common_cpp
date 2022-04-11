#ifndef SRC_MATRIX_MATRIXTYPES_H_
#define SRC_MATRIX_MATRIXTYPES_H_

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>
#include <Eigen/Sparse>

#if defined INCLUDE_NUMBER_THEORY_BOOST_GMP_INT ||                             \
    defined INCLUDE_NUMBER_THEORY_BOOST_CPP_INT
#include <boost/multiprecision/eigen.hpp>
#endif

template <typename T> using MyVector = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template <typename T>
using MyMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template <typename T>
using MySparseMatrix = Eigen::SparseMatrix<T, Eigen::ColMajor>;

#endif // SRC_MATRIX_MATRIXTYPES_H_
