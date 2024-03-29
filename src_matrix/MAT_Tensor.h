// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_TENSOR_H_
#define SRC_MATRIX_MAT_TENSOR_H_

// clang-format off
#include "MAT_Matrix.h"
#include <unsupported/Eigen/CXX11/Tensor>
// clang-format on

template <typename T>
MyMatrix<T> DimensionExtraction(Eigen::Tensor<T, 3> const &eT,
                                size_t const &iDim, int const &eDim) {
  int n1 = eT.dimension(0);
  int n2 = eT.dimension(1);
  int n3 = eT.dimension(2);
  if (iDim == 0) {
    MyMatrix<T> eMat(n2, n3);
    for (int i2 = 0; i2 < n2; i2++)
      for (int i3 = 0; i3 < n3; i3++)
        eMat(i2, i3) = eT(eDim, i2, i3);
    return eMat;
  }
  if (iDim == 1) {
    MyMatrix<T> eMat(n1, n3);
    for (int i1 = 0; i1 < n1; i1++)
      for (int i3 = 0; i3 < n3; i3++)
        eMat(i1, i3) = eT(i1, eDim, i3);
    return eMat;
  }
  if (iDim == 2) {
    MyMatrix<T> eMat(n1, n2);
    for (int i1 = 0; i1 < n1; i1++)
      for (int i2 = 0; i2 < n2; i2++)
        eMat(i1, i2) = eT(i1, i2, eDim);
    return eMat;
  }
  std::cerr << "Wrong input in ThreeDimArray\n";
  std::cerr << "iDim=" << iDim << "\n";
  std::cerr << "Allowed values: 0, 1, 2\n";
  throw TerminalException{1};
}

template <typename T, typename Foper>
T operationFromCoeff(Eigen::Tensor<T, 3> const &X, Foper f_oper) {
  auto LDim = X.dimensions();
  int a = LDim[0];
  int b = LDim[1];
  int c = LDim[2];
  T eRes = X(0, 0, 0);
  for (int i = 0; i < a; i++)
    for (int j = 0; j < b; j++)
      for (int k = 0; k < c; k++) {
        T eVal = X(i, j, k);
        f_oper(eRes, eVal);
      }
  return eRes;
}

template <typename T> T maxCoeff(Eigen::Tensor<T, 3> const &X) {
  auto f_max = [](T const &val1, T &val2) -> void {
    if (val1 > val2)
      val2 = val1;
  };
  return operationFromCoeff(X, f_max);
}

template <typename T> T minCoeff(Eigen::Tensor<T, 3> const &X) {
  auto f_max = [](T const &val1, T &val2) -> void {
    if (val1 < val2)
      val2 = val1;
  };
  return operationFromCoeff(X, f_max);
}

template <typename T>
Eigen::Tensor<T, 3> ZeroTensor3(int const &a, int const &b, int const &c) {
  Eigen::Tensor<T, 3> TheTens(a, b, c);
  for (int i = 0; i < a; i++)
    for (int j = 0; j < b; j++)
      for (int k = 0; k < c; k++)
        TheTens(i, j, k) = 0;
  return TheTens;
}

// clang-format off
#endif  // SRC_MATRIX_MAT_TENSOR_H_
// clang-format on
