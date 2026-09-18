// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// The Smith invariant factors of a matrix given in the sparse format of
// ReadSparseMatrix:
//
//     nbRow nbCol nnz
//     iRow iCol value
//     ...
//
// Same computation and same output as SmithNormalFormInvariant, the input
// representation being the only difference. It is the form to use for the
// boundary matrices of a chain complex, which are far too sparse to be
// worth holding densely: the dense pipeline costs time proportional to
// the number of cells, which for those matrices is three to four orders
// of magnitude more than the number of entries.
//
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheorySafeInt.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "MAT_MatrixIntSparse.h"
// clang-format on

template <typename T> void process(std::string const &FileI) {
  std::ifstream is(FileI);
  if (!is) {
    std::cerr << "Failed to open " << FileI << "\n";
    throw TerminalException{1};
  }
  MySparseMatrix<T> M = ReadSparseMatrix<T>(is);
  std::cerr << "|M|=" << M.rows() << " / " << M.cols()
            << " nnz=" << M.nonZeros() << "\n";
  MyVector<T> VectInv = SmithNormalFormInvariant_sparse(M);
  //
  std::map<T, size_t> MultInv;
  int len = VectInv.size();
  for (int u = 0; u < len; u++) {
    T val = VectInv(u);
    MultInv[val] += 1;
  }
  std::cerr << "MultInv =";
  for (auto &kv : MultInv) {
    std::cerr << " [" << kv.first << "," << kv.second << "]";
  }
  std::cerr << "\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "SmithNormalFormInvariant_sparse [arith] [inputMat]\n";
      std::cerr << "---\n";
      std::cerr << "arith: mpz_class, safe_integer, boost_cpp_int";
      std::cerr << ", flint_integer (with flint support)\n";
      std::cerr << "inputMat: the matrix in the sparse format\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    auto f = [&]() -> void {
      if (arith == "mpz_class")
        return process<mpz_class>(FileI);
      if (arith == "safe_integer")
        return process<SafeInt64>(FileI);
      if (arith == "boost_cpp_int")
        return process<boost::multiprecision::cpp_int>(FileI);
#ifdef ENABLE_FLINT_SUPPORT
      if (arith == "flint_integer")
        return process<fmpz_class>(FileI);
#endif
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
