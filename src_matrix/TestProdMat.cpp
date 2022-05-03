#include "MAT_Matrix.h"
#include "NumberTheory.h"
#include "rational.h"
int main(int argc, char *argv[]) {
  //  using T=mpq_class;
  using T = Rational<long>;
  try {
    if (argc != 3) {
      fprintf(stderr, "TestProdMat [M1] [M2]\n");
      return -1;
    }
    // reading the matrix
    std::ifstream is1(argv[1]);
    std::ifstream is2(argv[2]);
    MyMatrix<T> M1 = ReadMatrix<T>(is1);
    MyMatrix<T> M2 = ReadMatrix<T>(is2);
    std::cerr << "M1=\n";
    WriteMatrix(std::cerr, M1);
    std::cerr << "M2=\n";
    WriteMatrix(std::cerr, M2);
    //
    MyMatrix<T> eProd = M1 * M2;
    std::cerr << "eProd=\n";
    WriteMatrix(std::cerr, eProd);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}