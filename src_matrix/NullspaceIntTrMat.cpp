// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixInt.h"
// clang-format on

template<typename T>
void process(std::string const& FileI, std::ostream & os) {
  MyMatrix<T> TheMat = ReadMatrixFile<T>(FileI);
  // computing the kernel
  MyMatrix<T> KerInt = NullspaceIntTrMat(TheMat);
  os << "return ";
  WriteMatrixGAP(os, KerInt);
  os << ";\n";
}


int main(int argc, char *argv[]) {
  try {
    if (argc != 4 && argc != 3) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "NullspaceIntTrMat [arith] [inputMat] [output]\n");
      fprintf(stderr, "or\n");
      fprintf(stderr, "NullspaceInttrMat [arith] [inputMat]\n");
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::ifstream INmat(argv[1]);

    auto Prt = [&](std::ostream &os) -> void {
      if (arith == "mpz_class")
        process<mpz_class>(FileI, os);
      if (arith == "mpq_class")
        process<mpq_class>(FileI, os);
      if (arith == "safe_integer")
        process<SafeInt64>(FileI, os);
      if (arith == "safe_rational")
        process<Rational<SafeInt64>>(FileI, os);
      std::cerr << "Failed to find a matching type\n";
      throw TerminalException{1};
    };
    //
    if (argc == 4) {
      std::string FileO = argv[3];
      std::ofstream os(argv[2]);
      Prt(os);
    } else {
      Prt(std::cerr);
    }
    std::cerr << "Normal termination of NullspaceIntTrMat\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong for NullspaceIntTrMat\n";
    exit(e.eVal);
  }
}
