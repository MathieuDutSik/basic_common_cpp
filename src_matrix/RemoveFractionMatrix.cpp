#include "MAT_MatrixInt.h"
#include "NumberTheory.h"
int main(int argc, char *argv[]) {
  using T = mpq_class;
  try {
    if (argc != 2) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "RemoveFractionMatrix [inputMat]\n");
      return -1;
    }
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> M = ReadMatrix<T>(INmat);
    //
    std::cerr << "M=\n";
    WriteMatrix(std::cerr, M);
    //
    FractionMatrix<T> eRec = RemoveFractionMatrixPlusCoeff(M);
    //
    std::cerr << "Mult=" << eRec.TheMult << "\n";
    std::cerr << "TheMat=\n";
    WriteMatrix(std::cerr, eRec.TheMat);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
