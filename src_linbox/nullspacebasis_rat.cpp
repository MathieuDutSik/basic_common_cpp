// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

#include "linbox/algorithms/gauss.h"
#include "linbox/matrix/dense-matrix.h"
#include <givaro/gfq.h>
#include <iostream>

using namespace LinBox;

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Usage to get a null space basis over Q:  "
                 "<matrix-file-in-SMS-format>"
              << std::endl;
    return -1;
  }

  std::ifstream input(argv[1]);
  if (!input) {
    std::cerr << "Error opening matrix file " << argv[1] << std::endl;
    return -1;
  }

  using Rats = Givaro::QField<Givaro::Rational>;
  Rats QQ;
  SparseMatrix<Rats, SparseMatrixFormat::SparseSeq> B(QQ);
  B.read(input);
  std::cout << "B is " << B.rowdim() << " by " << B.coldim() << std::endl;

  DenseMatrix<Rats> NullSpace(QQ, B.coldim(), B.coldim());
  GaussDomain<Rats> GD(QQ);

  GD.nullspacebasisin(NullSpace, B);

  NullSpace.write(std::cerr << "X:=", Tag::FileFormat::Maple)
      << ';' << std::endl;

  std::cerr << "NullsSpace dimensions:" << NullSpace.rowdim() << 'x'
            << NullSpace.coldim() << std::endl;

  return 0;
}
