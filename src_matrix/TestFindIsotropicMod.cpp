// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
#include "MAT_MatrixMod.h"
// clang-format on

template <typename T>
void find_isotropic_kernel(std::string const &FileI) {
  std::ifstream is(FileI);
  MyMatrix<T> M = ReadMatrix<T>(is);
  T TheMod;
  is >> TheMod;
  //
  std::optional<MyVector<T>> opt = FindIsotropicVectorMod(M, TheMod);
  if (opt) {
    MyVector<T> const& V = *opt;
    std::cerr << "Found mod isotropic vector V=" << StringVectorGAP(V) << "\n";
  } else {
    std::cerr << "Failed to find an isotropic vector\n";
  }

}

void find_isotropic(std::string const &arithmetic,
                    std::string const &FileI) {
  if (arithmetic == "safe_integer") {
    using T = SafeInt64;
    return find_isotropic_kernel<T>(FileI);
  }
  if (arithmetic == "integer") {
    using T = mpz_class;
    return find_isotropic_kernel<T>(FileI);
  }
  if (arithmetic == "rational") {
    using T = mpq_class;
    return find_isotropic_kernel<T>(FileI);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "SolutionIntMat [arithmetic] [inputMat]\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string FileI = argv[2];
    find_isotropic(arithmetic, FileI);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
