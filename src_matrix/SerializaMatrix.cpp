#include "MAT_Matrix.h"
#include "MAT_MatrixSerialization.h"
#include "NumberTheory.h"


int main() {
  //  using T = mpz_class;
  using T = int;
  size_t n_row = 10;
  size_t n_col = 20;
  MyMatrix<T> M1(n_row, n_col);
  for (size_t i_row=0; i_row<n_row; i_row++) {
    for (size_t i_col=0; i_col<n_col; i_col++) {
      T val = rand() % 20;
      M1(i_row, i_col) = val;
    }
  }
  std::string filename = "/tmp/MAT_filename.boost_archive";

  //
  // Writing the data 
  //

  // save data to archive
  {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);
    // write class instance to archive
    oa << M1;
    // archive and stream closed when destructors are called
  }

  MyMatrix<T> M2;
  // load data from archive
  {
    // create and open an archive for input
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);
    // read class state from archive
    ia >> M2;
    // archive and stream closed when destructors are called
  }

  // Checking the consistency of the data exchange
  size_t hash1 = std::hash<MyMatrix<T>>()(M1);
  size_t hash2 = std::hash<MyMatrix<T>>()(M2);
  std::cerr << "hash1=" << hash1 << " hash2=" << hash2 << "\n";
}
