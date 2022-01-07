//#include "NumberTheoryBoost.h"
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"



int main(int argc, char *argv[])
{
  //  using T=mpq_class;
  using T=mpz_class;
  //  using T=boost::multiprecision::cpp_int;
  //  using T=int;
  //  using T=long;
  try {
    if (argc != 3) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "RowHermiteNormalForm [inputMat] [output]\n");
      return -1;
    }
    // reading the matrix
    std::ifstream INmat(argv[1]);
    MyMatrix<T> TheMat=ReadMatrix<T>(INmat);
    // computing the kernel
    std::pair<MyMatrix<T>, MyMatrix<T>> ePair = ComputeRowHermiteNormalForm(TheMat);
    //    MyMatrix<T> TheKer=NullspaceMat(TheMat);
    //
    std::ofstream os(argv[2]);
    os << "return [";
    WriteMatrixGAP(os, ePair.first);
    os << ",\n";
    WriteMatrixGAP(os, ePair.second);
    os << "];\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
