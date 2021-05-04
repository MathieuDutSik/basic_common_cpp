#ifndef INCLUDE_MATRIX_TYPES_H
#define INCLUDE_MATRIX_TYPES_H




#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include <Eigen/Eigenvalues>

template <typename T>
using MyVector = Eigen::Matrix<T,Eigen::Dynamic,1>;

template <typename T>
using MyMatrix = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>;

template <typename T>
using MySparseMatrix = Eigen::SparseMatrix<T,Eigen::ColMajor>;

#endif
