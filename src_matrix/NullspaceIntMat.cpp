#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
int main(int argc, char *argv[])
{
  try {
    if (argc != 3 && argc != 2) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "NullspaceIntMat [inputMat] [output]\n");
      fprintf(stderr, "or\n");
      fprintf(stderr, "NullspaceIntMat [inputMat]\n");
      return -1;
    }
    //using T=mpz_class;
    using T=mpq_class;
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> TheMat=ReadMatrix<T>(INmat);
    // computing the kernel
    MyMatrix<T> KerInt = NullspaceIntMat(TheMat);
    MyMatrix<T> TheKer = NullspaceMat(TheMat);
    if (KerInt.cols() != TheKer.cols() || KerInt.rows() != TheKer.rows()) {
      std::cerr << "Difference in dimension between NullspaceIntMat and NullspaceMat\n";
      throw TerminalException{1};
    }
    //
    auto Prt=[&](std::ostream & os) -> void {
      os << "return ";
      WriteMatrixGAP(os, KerInt);
      os << ";\n";
    };
    //
    if (argc == 3) {
      std::ofstream os(argv[2]);
      Prt(os);
    } else {
      Prt(std::cerr);
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
