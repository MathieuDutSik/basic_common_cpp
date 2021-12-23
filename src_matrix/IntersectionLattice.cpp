#include "MAT_MatrixInt.h"
#include "NumberTheory.h"
int main(int argc, char *argv[])
{
  //  using T=mpz_class;
  using T=mpq_class;
  try {
    int n=2;
    auto get_fullrank=[&]() -> MyMatrix<T> {
      MyMatrix<T> M(n,n);
      while(true) {
        for (int i=0; i<n; i++)
          for (int j=0; j<n; j++)
            M(i,j) = rand() % 7;
        if (RankMat(M) == n)
          return M;
      }
    };
    for (int iter=0; iter<10; iter++) {
      std::cerr << "iter=" << iter << "\n";
      MyMatrix<T> M1 = get_fullrank();
      MyMatrix<T> M2 = get_fullrank();
      MyMatrix<T> M1_i_M2 = IntersectionLattice(M1, M2);
      std::cerr << "   M1=\n";
      WriteMatrix(std::cerr, M1);
      std::cerr << "   M2=\n";
      WriteMatrix(std::cerr, M2);
      std::cerr << "   M1_i_M2=\n";
      WriteMatrix(std::cerr, M1_i_M2);
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
