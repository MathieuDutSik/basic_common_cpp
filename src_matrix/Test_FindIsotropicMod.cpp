// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
#include "MAT_MatrixMod.h"
// clang-format on

template <typename T> void find_isotropic_kernel(std::string const &FileI) {
  std::ifstream is(FileI);
  MyMatrix<T> M = ReadMatrix<T>(is);
  T TheMod;
  is >> TheMod;
  //
  std::optional<MyVector<T>> opt = FindIsotropicVectorMod(M, TheMod);
  if (opt) {
    MyVector<T> const &V = *opt;
    std::cerr << "Found mod isotropic vector V=" << StringVectorGAP(V) << "\n";
  } else {
    std::cerr << "Failed to find an isotropic vector\n";
  }
}

void find_isotropic(std::string const &arithmetic, std::string const &FileI) {
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
  if (arithmetic == "safe_rational") {
    using T = Rational<SafeInt64>;
    return find_isotropic_kernel<T>(FileI);
  }
  if (arithmetic == "boost_cpp_int") {
    using T = boost::multiprecision::cpp_int;
    return find_isotropic_kernel<T>(FileI);
  }
  if (arithmetic == "boost_cpp_rational") {
    using T = boost::multiprecision::cpp_rational;
    return find_isotropic_kernel<T>(FileI);
  }
  if (arithmetic == "boost_mpz_int") {
    using T = boost::multiprecision::mpz_int;
    return find_isotropic_kernel<T>(FileI);
  }
  if (arithmetic == "boost_mpq_rational") {
    using T = boost::multiprecision::mpq_rational;
    return find_isotropic_kernel<T>(FileI);
  }
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "Test_FindIsotropicMod [arithmetic] [inputMat]\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string FileI = argv[2];
    find_isotropic(arithmetic, FileI);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
