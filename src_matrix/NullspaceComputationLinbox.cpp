#include "NumberTheory.h"
#include "MatrixLinbox.h"
int main(int argc, char *argv[])
{
  using T=mpq_class;
  try {
    if (argc != 3) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "NullspaceComputationLinbox [inputMat] [output]\n");
      return -1;
    }
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> TheMat = ReadMatrix<T>(INmat);
    MyMatrix<T> TheKer = NullspaceTrMat_linbox(TheMat);
    //
    std::ofstream os(argv[2]);
    os << "return ";
    WriteMatrixGAP(os, TheKer);
    os << ";\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
