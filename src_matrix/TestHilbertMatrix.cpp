#include "NumberTheory.h"
#include "rational.h"
#include "MAT_Matrix.h"
int main(int argc, char *argv[])
{
  //  using T=mpq_class;
  using T=Rational<long>;
  try {
    if (argc != 2) {
      fprintf(stderr, "TestHilbertMatrix is used as\n");
      fprintf(stderr, "TestHilbertMatrix [n]\n");
      return -1;
    }
    // reading the matrix
    int n;
    sscanf(argv[1], "%d", &n);
    //
    MyMatrix<T> eMat(n,n);
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        T val = i + j + 1;
        eMat(i,j) = 1 / val;
      }
    }
    std::cerr << "eMat=\n";
    WriteMatrix(std::cerr, eMat);

    
    MyMatrix<T> eInv = Inverse(eMat);
    std::cerr << "eInv=\n";
    WriteMatrix(std::cerr, eInv);
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
