#include "NumberTheory.h"
#include "MAT_Matrix.h"
int main(int argc, char *argv[])
{
  using T=mpq_class;
  //  using T=mpz_class;
  try {
    if (argc != 2) {
      fprintf(stderr, "TestHilbertMatrix is used as\n");
      fprintf(stderr, "TestHilbertMatrix [n]\n");
      return -1;
    }
    // reading the matrix
    int n;
    sscanf(argv[1], "%d", &n);
    int siz = 2;
    MyMatrix<T> eMat(n,m);
    for (int i=0; i<n; i++) {
      for (int j=0; j<m; j++) {
        T val = i + j + 1;
        eMat(i,j) = 1 / val;
      }
    }
    MyMatrix<T> eInv = Inverse(eMat);
    std::cerr << "eInv=\n";
    WriteMatrix(std::cerr, eInv);
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
