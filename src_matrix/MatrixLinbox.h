#ifndef INCLUDE_MATRIX_LINBOX_H
#define INCLUDE_MATRIX_LINBOX_H

#include "MAT_Matrix.h"
#include "NumberGivaro.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/algorithms/gauss.h"



MyMatrix<mpq_class> NullspaceTrMat_linbox(MyMatrix<mpq_class> const& M)
{
  using Rats=Givaro::QField<Givaro::Rational>;
  Rats QQ;
  size_t n_rows = M.rows();
  size_t n_cols = M.cols();

  LinBox::SparseMatrix<Rats> B(QQ);
  B.resize(n_rows, n_cols);
  for (size_t i_row=0; i_row<n_rows; i_row++)
    for (size_t i_col=0; i_col<n_cols; i_col++) {
      Givaro::Rational val = GetGivaroRational(M(i_row, i_col));
      B.appendEntry(i_row, i_col, val);
    }
  // The NullSpace matrix. Maybe we can simplify to smaller matrix
  // but so far no success.
  LinBox::DenseMatrix<Rats> NullSpace(QQ, n_cols, n_cols);
  LinBox::GaussDomain<Rats> GD(QQ);

  GD.nullspacebasisin(NullSpace, B);
  Givaro::Rational zero(0/1);
  auto iszero=[&](size_t i_row) -> bool {
    for (size_t i_col=0; i_col<n_cols; i_col++) {
      Givaro::Rational val = NullSpace.getEntry(i_col, i_row);
      if (val != zero)
        return false;
    }
    return true;
  };
  size_t rank = 0;
  while(true) {
    if (iszero(rank))
      break;
    rank++;
  }
  //
  MyMatrix<mpq_class> Ker(rank, n_cols);
  for (size_t i_row=0; i_row<rank; i_row++) {
    for (size_t i_col=0; i_col<n_cols; i_col++) {
      Givaro::Rational val = NullSpace.getEntry(i_col, i_row);
      Ker(i_row, i_col) = ConvertGivaroRational(val);
    }
  }
  return Ker;
}
















#endif
