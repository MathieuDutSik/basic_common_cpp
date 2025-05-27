// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "MAT_MatrixInt.h"
// clang-format on

template <typename T> void process(std::string const &FileI) {
  MyMatrix<T> M = ReadMatrixFile<T>(FileI);
  std::cerr << "|M|=" << M.rows() << " / " << M.cols() << "\n";
  // computing the Smith normal form
  MyVector<T> VectInv = SmithNormalFormInvariant(M);
  //
  //  std::cerr << "VectInv=\n";
  //  WriteVector(std::cerr, VectInv);
  //
  std::map<T, size_t> MultInv;
  int len = VectInv.size();
  for (int u=0; u<len; u++) {
    T val = VectInv(u);
    MultInv[val] += 1;
  }
  std::cerr << "MultInv =";
  for (auto & kv: MultInv) {
    std::cerr << " [" << kv.first << "," << kv.second << "]";
  }
  std::cerr << "\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "SmithNormalForm [arith] [inputMat]\n";
      std::cerr << "---\n";
      std::cerr << "arith: mpz_class, mpq_class, safe_integer, safe_rational\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    auto f = [&]() -> void {
      if (arith == "mpz_class")
        return process<mpz_class>(FileI);
      if (arith == "mpq_class")
        return process<mpq_class>(FileI);
      if (arith == "safe_integer")
        return process<SafeInt64>(FileI);
      if (arith == "safe_rational")
        return process<Rational<SafeInt64>>(FileI);
      std::cerr << "Failed to find a matching type\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
