#include "MAT_MatrixInt.h"
#include "NumberTheory.h"
int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "NullspaceComputation [inputMat] [output]\n");
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
    std::ofstream os(argv[2]);
    os << "return ";
    WriteMatrixGAP(os, KerInt);
    os << ";\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
