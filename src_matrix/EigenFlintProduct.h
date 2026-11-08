// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_EIGENFLINTPRODUCT_H_
#define SRC_MATRIX_EIGENFLINTPRODUCT_H_

// The Eigen matrix products of fmpz_class / fmpq_class computed by
// fmpz_mat_mul / fmpq_mat_mul. Eigen evaluates every matrix-matrix product
// through Eigen::internal::general_matrix_matrix_product and every
// matrix-vector product through general_matrix_vector_product, two class
// templates on the scalar type that receive the operands as pointers and
// strides. Their partial specializations below, constrained to the flint
// types, replace the generic scalar-by-scalar kernel, as the BLAS backend of
// Eigen does for float and double. So a * b * c, a * b.transpose(),
// a.block(...) * v and so on reach flint with no change at the call site.
//
// No copy: fmpz_class is one fmpz and fmpq_class one fmpq, and a flint
// matrix is a pointer to row-major entries with a row stride, so a flint
// matrix can be a view on the Eigen storage. A column-major Eigen matrix
// viewed that way is the transpose, which is taken into account by
// computing res^T = rhs^T * lhs^T. An operand in the other storage order
// (a transposed view, a row-major matrix) is copied once, transposed.
//
// The products of rows + cols + depth < 20 go there too through
// eigen_gemm_coeff_threshold (EigenProductThreshold.h), except with Eigen
// 3.3 which does not have that macro.
//
// This header is included by MatrixTypes.h with ENABLE_FLINT_SUPPORT.

// clang-format off
#include "NumberTheoryFlint.h"
#include <flint/fmpz_mat.h>
#include <flint/fmpq_mat.h>
// clang-format on

static_assert(sizeof(fmpz_class) == sizeof(fmpz),
              "fmpz_class must be a single fmpz for the flint matrix views");
static_assert(sizeof(fmpq_class) == sizeof(fmpq),
              "fmpq_class must be a single fmpq for the flint matrix views");
static_assert(requires(fmpz_mat_struct m) { m.stride; },
              "the flint matrix views need flint >= 3.1 (fmpz_mat_struct "
              "with a row stride)");

template <> struct eigen_gemm_coeff_threshold<fmpz_class> {
  static constexpr int value = 0;
};

template <> struct eigen_gemm_coeff_threshold<fmpq_class> {
  static constexpr int value = 0;
};

template <typename T> struct flint_mat_ops;

template <> struct flint_mat_ops<fmpz_class> {
  using entry = fmpz;
  using mat = fmpz_mat_struct;
  static void init(mat *m, slong r, slong c) { fmpz_mat_init(m, r, c); }
  static void clear(mat *m) { fmpz_mat_clear(m); }
  static void set_entry(mat *m, slong i, slong j, fmpz_class const &x) {
    fmpz_set(fmpz_mat_entry(m, i, j), x.get_fmpz_t());
  }
  static void transpose(mat *b, mat const *a) { fmpz_mat_transpose(b, a); }
  static void mul(mat *c, mat const *a, mat const *b) { fmpz_mat_mul(c, a, b); }
  static void add(mat *c, mat const *a, mat const *b) { fmpz_mat_add(c, a, b); }
  static void addmul_scalar(mat *c, mat const *a, fmpz_class const &x) {
    fmpz_mat_scalar_addmul_fmpz(c, a, x.get_fmpz_t());
  }
  static bool is_one(fmpz_class const &x) { return fmpz_is_one(x.get_fmpz_t()); }
  static bool is_zero(mat const *m) { return fmpz_mat_is_zero(m); }
};

template <> struct flint_mat_ops<fmpq_class> {
  using entry = fmpq;
  using mat = fmpq_mat_struct;
  static void init(mat *m, slong r, slong c) { fmpq_mat_init(m, r, c); }
  static void clear(mat *m) { fmpq_mat_clear(m); }
  static void set_entry(mat *m, slong i, slong j, fmpq_class const &x) {
    fmpq_set(fmpq_mat_entry(m, i, j), x.get_fmpq_t());
  }
  static void transpose(mat *b, mat const *a) { fmpq_mat_transpose(b, a); }
  static void mul(mat *c, mat const *a, mat const *b) { fmpq_mat_mul(c, a, b); }
  static void add(mat *c, mat const *a, mat const *b) { fmpq_mat_add(c, a, b); }
  static void addmul_scalar(mat *c, mat const *a, fmpq_class const &x) {
    mat t;
    fmpq_mat_init(&t, a->r, a->c);
    fmpq_mat_scalar_mul_fmpq(&t, a, x.get_fmpq_t());
    fmpq_mat_add(c, c, &t);
    fmpq_mat_clear(&t);
  }
  static bool is_one(fmpq_class const &x) { return fmpq_is_one(x.get_fmpq_t()); }
  static bool is_zero(mat const *m) { return fmpq_mat_is_zero(m); }
};

// The r x c flint matrix of row-major entries at p, rows `stride` apart.
// The entries stay owned by the Eigen matrix.
template <typename T>
typename flint_mat_ops<T>::mat flint_view(T const *p, slong r, slong c,
                                          slong stride) {
  typename flint_mat_ops<T>::mat m;
  m.entries =
      reinterpret_cast<typename flint_mat_ops<T>::entry *>(const_cast<T *>(p));
  m.r = r;
  m.c = c;
  m.stride = stride;
  return m;
}

// res += alpha * a * b. Eigen zeroes the destination before evaluating a
// product into it, and then the product is written into res directly.
template <typename T>
void flint_mul_add(typename flint_mat_ops<T>::mat *res,
                   typename flint_mat_ops<T>::mat const *a,
                   typename flint_mat_ops<T>::mat const *b, T const &alpha) {
  using Ops = flint_mat_ops<T>;
  if (Ops::is_one(alpha) && Ops::is_zero(res)) {
    Ops::mul(res, a, b);
    return;
  }
  typename Ops::mat prod;
  Ops::init(&prod, res->r, res->c);
  Ops::mul(&prod, a, b);
  if (Ops::is_one(alpha))
    Ops::add(res, res, &prod);
  else
    Ops::addmul_scalar(res, &prod, alpha);
  Ops::clear(&prod);
}

// The transpose of the r x c matrix stored at p in the order `Order`, as a
// flint matrix: a view for a column-major matrix, a copy otherwise.
template <typename T, int Order> struct flint_transposed_operand {
  using Ops = flint_mat_ops<T>;
  typename Ops::mat m;
  flint_transposed_operand(T const *p, slong r, slong c, slong stride) {
    if constexpr (Order == Eigen::ColMajor) {
      m = flint_view(p, c, r, stride);
    } else {
      typename Ops::mat v = flint_view(p, r, c, stride);
      Ops::init(&m, c, r);
      Ops::transpose(&m, &v);
    }
  }
  ~flint_transposed_operand() {
    if constexpr (Order != Eigen::ColMajor)
      Ops::clear(&m);
  }
  flint_transposed_operand(flint_transposed_operand const &) = delete;
  flint_transposed_operand &
  operator=(flint_transposed_operand const &) = delete;
};

// res (rows x cols, column-major) += alpha * lhs (rows x depth) * rhs
// (depth x cols), computed as res^T += alpha * rhs^T * lhs^T.
template <typename T, int LhsOrder, int RhsOrder>
void flint_gemm(slong rows, slong cols, slong depth, T const *lhs,
                slong lhsStride, T const *rhs, slong rhsStride, T *res,
                slong resStride, T const &alpha) {
  flint_transposed_operand<T, LhsOrder> lhsT(lhs, rows, depth, lhsStride);
  flint_transposed_operand<T, RhsOrder> rhsT(rhs, depth, cols, rhsStride);
  typename flint_mat_ops<T>::mat resT = flint_view(res, cols, rows, resStride);
  flint_mul_add<T>(&resT, &rhsT.m, &lhsT.m, alpha);
}

// res (rows entries, resIncr apart) += alpha * lhs (rows x cols) * rhs (cols
// entries, rhsIncr apart), as a flint product with a one-row or one-column
// matrix. A column-major lhs is viewed as its transpose and the product
// computed as res^T = rhs^T * lhs^T, which needs rhs as a contiguous row.
template <typename T, int LhsOrder>
void flint_gemv(slong rows, slong cols, T const *lhs, slong lhsStride,
                T const *rhs, slong rhsIncr, T *res, slong resIncr,
                T const &alpha) {
  using Ops = flint_mat_ops<T>;
  using mat = typename Ops::mat;
  if constexpr (LhsOrder == Eigen::ColMajor) {
    // Eigen passes a contiguous result in that case.
    eigen_internal_assert(resIncr == 1);
    mat lhsT = flint_view(lhs, cols, rows, lhsStride);
    mat r = flint_view(res, 1, rows, rows);
    if (rhsIncr == 1) {
      mat x = flint_view(rhs, 1, cols, cols);
      flint_mul_add<T>(&r, &x, &lhsT, alpha);
    } else {
      mat x;
      Ops::init(&x, 1, cols);
      for (slong k = 0; k < cols; k++)
        Ops::set_entry(&x, 0, k, rhs[k * rhsIncr]);
      flint_mul_add<T>(&r, &x, &lhsT, alpha);
      Ops::clear(&x);
    }
  } else {
    mat a = flint_view(lhs, rows, cols, lhsStride);
    mat x = flint_view(rhs, cols, 1, rhsIncr);
    mat r = flint_view(res, rows, 1, resIncr);
    flint_mul_add<T>(&r, &a, &x, alpha);
  }
}

namespace Eigen::internal {

// Eigen calls it with a column-major result after swapping the operands of
// a row-major one, and with an inner stride 1 for a plain result.
template <typename Index, typename T, int LhsOrder, int RhsOrder>
  requires is_flint_type<T>::value
struct general_matrix_matrix_product<Index, T, LhsOrder, false, T, RhsOrder,
                                     false, ColMajor, 1> {
  typedef gebp_traits<T, T> Traits;
  static void run(Index rows, Index cols, Index depth, const T *lhs,
                  Index lhsStride, const T *rhs, Index rhsStride, T *res,
                  Index /*resIncr*/, Index resStride, T alpha,
                  level3_blocking<T, T> & /*blocking*/,
                  GemmParallelInfo<Index> * /*info*/ = 0) {
    flint_gemm<T, LhsOrder, RhsOrder>(rows, cols, depth, lhs, lhsStride, rhs,
                                      rhsStride, res, resStride, alpha);
  }
};

// One specialization per storage order, as Eigen has: a single one generic
// in the order would be ambiguous with them. The mappers carry the pointer
// and the stride (for the right side, the increment between its entries).
template <typename Index, typename T, typename LhsMapper, typename RhsMapper,
          int Version>
  requires is_flint_type<T>::value
struct general_matrix_vector_product<Index, T, LhsMapper, ColMajor, false, T,
                                     RhsMapper, false, Version> {
  static void run(Index rows, Index cols, const LhsMapper &lhs,
                  const RhsMapper &rhs, T *res, Index resIncr, T alpha) {
    flint_gemv<T, ColMajor>(rows, cols, lhs.data(), lhs.stride(), rhs.data(),
                            rhs.stride(), res, resIncr, alpha);
  }
};

template <typename Index, typename T, typename LhsMapper, typename RhsMapper,
          int Version>
  requires is_flint_type<T>::value
struct general_matrix_vector_product<Index, T, LhsMapper, RowMajor, false, T,
                                     RhsMapper, false, Version> {
  static void run(Index rows, Index cols, const LhsMapper &lhs,
                  const RhsMapper &rhs, T *res, Index resIncr, T alpha) {
    flint_gemv<T, RowMajor>(rows, cols, lhs.data(), lhs.stride(), rhs.data(),
                            rhs.stride(), res, resIncr, alpha);
  }
};

}  // namespace Eigen::internal

// clang-format off
#endif  // SRC_MATRIX_EIGENFLINTPRODUCT_H_
// clang-format on
