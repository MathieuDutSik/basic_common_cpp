#include "MAT_MatrixInt.h"
#include "NumberTheory.h"
int main(int argc, char *argv[]) {
  try {
    if (argc != 2) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "SmithNormalForm [inputMat]\n");
      return -1;
    }
    // using T=mpz_class;
    using T = mpq_class;
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> M = ReadMatrix<T>(INmat);
    // computing the kernel
    std::pair<MyMatrix<T>, MyMatrix<T>> pair = SmithNormalForm(M);
    std::cerr << "ROW=\n";
    WriteMatrix(std::cerr, pair.first);
    std::cerr << "COL=\n";
    WriteMatrix(std::cerr, pair.second);
    //
    MyMatrix<T> RedMat = pair.first * M * pair.second;
    std::cerr << "RedMat=\n";
    WriteMatrix(std::cerr, RedMat);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
