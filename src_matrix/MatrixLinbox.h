#ifndef INCLUDE_MATRIX_LINBOX_H
#define INCLUDE_MATRIX_LINBOX_H

#include <MAT_Matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include <givaro/gfq.h>
#include "linbox/algorithms/gauss.h"



MyMatrix<mpq_class> NullspaceTrMat_linbox(MyMatrix<mpq_class> const& M, size_t expected_rank)
{
  using Rats=Givaro::QField<Givaro::Rational>;
  Rats QQ;
  size_t n_rows = M.rows();
  size_t n_cols = M.cols();

  LinBox::DenseMatrix<Rats> B(QQ, n_rows, n_cols);
  for (size_t i_row=0; i_row<n_rows; i_row++)
    for (size_t i_col=0; i_col<n_cols; i_col++) {
      Givaro::Rational val = GetGivaroRational(M(i_row, i_col));
      B.setEntry(i,j, val);
    }
  DenseMatrix<Rats> NullSpace(QQ, n_cols, expected_rank);
  GaussDomain<Rats> GD(QQ);

  GD.nullspacebasisin(NullSpace, B);
  auto iszero=[&](size_t i_row) -> bool {
    for (size_t i_col=0; i_col<n_cols; i_col++) {
      Givaro::Rational val = NullSpace.getEntry(i_row, i_col);
      if (val != 0)
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
  MyMatrix<mpq_class> M(rank, n_cols);
  for (size_t i_row=0; i_row<rank; i_row++) {
    for (size_t i_col=0; i_col<n_cols; i_col++) {
      Givaro::Rational val = NullSpace.getEntry(i_row, i_col);
      M(i_row, i_col) = ConvertGivaroRational(val);
    }
  }
  return M;
}
















#endif
