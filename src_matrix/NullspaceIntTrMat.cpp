#include "MAT_MatrixInt.h"
#include "NumberTheory.h"
int main(int argc, char *argv[])
{
  try {
    if (argc != 3 && argc != 2) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "NullspaceIntTrMat [inputMat] [output]\n");
      fprintf(stderr, "or\n");
      fprintf(stderr, "NullspaceInttrMat [inputMat]\n");
      return -1;
    }
    using Tint=mpz_class;
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<Tint> TheMat=ReadMatrix<Tint>(INmat);
    // computing the kernel
    MyMatrix<Tint> KerInt = NullspaceIntTrMat(TheMat);
    //    MyMatrix<T> TheKer=NullspaceMat(TheMat);
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
