#ifndef INCLUDE_MATRIX_TYPES_H
#define INCLUDE_MATRIX_TYPES_H




#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>

#if defined INCLUDE_NUMBER_THEORY_BOOST_GMP_INT || defined INCLUDE_NUMBER_THEORY_BOOST_CPP_INT
#include <boost/multiprecision/eigen.hpp>
#endif


template <typename T>
using MyVector = Eigen::Matrix<T,Eigen::Dynamic,1>;

template <typename T>
using MyMatrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;

template <typename T>
using MySparseMatrix = Eigen::SparseMatrix<T,Eigen::ColMajor>;

#endif
