// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "MAT_Matrix.h"
// clang-format off

template <typename Tfull> void test_type(Tfull const& val1) {
  std::string filename = "/tmp/MAT_filename.boost_archive";

  //
  // Writing the data
  //

  // save data to archive
  {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    // write class instance to archive
    oa << val1;
    // archive and stream closed when destructors are called
  }

  Tfull val2;
  // load data from archive
  {
    // create and open an archive for input
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    // read class state from archive
    ia >> val2;
    // archive and stream closed when destructors are called
  }

  // Checking the consistency of the data exchange
  if (val1 != val2) {
    std::cerr << "Error in the serialization\n";
    throw TerminalException{1};
  }
}






template <typename T> void test_type_mymatrix_myvector() {
  size_t n_row = 10;
  size_t n_col = 20;
  MyMatrix<T> M(n_row, n_col);
  for (size_t i_row = 0; i_row < n_row; i_row++) {
    for (size_t i_col = 0; i_col < n_col; i_col++) {
      T val = UniversalScalarConversion<T,long>(random() % 20);
      M(i_row, i_col) = val;
    }
  }
  test_type(M);
  //
  MyVector<T> V(n_row);
  for (size_t i_row = 0; i_row < n_row; i_row++) {
    T val = UniversalScalarConversion<T,long>(random() % 20);
    V(i_row) = val;
  }
  test_type(V);
}

int main() {
  test_type_mymatrix_myvector<mpz_class>();
  test_type_mymatrix_myvector<mpq_class>();
  test_type_mymatrix_myvector<SafeInt64>();
  test_type_mymatrix_myvector<Rational<SafeInt64>>();
}
