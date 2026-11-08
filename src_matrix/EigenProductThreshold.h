// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_EIGENPRODUCTTHRESHOLD_H_
#define SRC_MATRIX_EIGENPRODUCTTHRESHOLD_H_

// Eigen evaluates a matrix product with rows + cols + depth below
// EIGEN_GEMM_TO_COEFFBASED_THRESHOLD (20 by default) coefficient by
// coefficient, without calling the GEMM kernel. The macro is expanded inside
// the product evaluator, where the scalar type is in scope as Scalar, so it
// can be made to depend on the scalar: the types whose GEMM kernel is
// replaced by a native one (the flint types, see EigenFlintProduct.h) set it
// to 0 so that every product of theirs reaches that kernel.
//
// The macro has to be defined before Eigen is included, so this header is
// included before every Eigen header. Eigen 3.3 does not have the macro
// (the threshold is the literal 20) and the small products stay
// coefficient-based there.

#if defined(EIGEN_CORE_H) || defined(EIGEN_CORE_MODULE_H)
#error "Eigen was included before EigenProductThreshold.h: include MatrixTypes.h first"
#endif

template <typename T> struct eigen_gemm_coeff_threshold {
  static constexpr int value = 20;
};

#define EIGEN_GEMM_TO_COEFFBASED_THRESHOLD                                     \
  (eigen_gemm_coeff_threshold<Scalar>::value)

// clang-format off
#endif  // SRC_MATRIX_EIGENPRODUCTTHRESHOLD_H_
// clang-format on
