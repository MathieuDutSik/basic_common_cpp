// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
#include "MAT_MatrixInt.h"
// clang-format on

template<typename T>
void compute_determinant_kernel(std::string const& eFile) {
  MyMatrix<T> A = ReadMatrixFile<T>(eFile);
  std::cerr << "Computing determinant of matrix A=\n";
  WriteMatrix(std::cerr, A);
  T TheDet_gauss = DeterminantMat(A);
  T TheDet_symm = DeterminantMatPermutation(A);
  std::cerr << "TheDet_gauss=" << TheDet_gauss << "\n";
  std::cerr << "TheDet_symm =" << TheDet_symm << "\n";
}

void compute_determinant(std::string const& arithmetic, std::string const& eFile) {
  if (arithmetic == "safe_integer") {
    using T = SafeInt64;
    return compute_determinant_kernel<T>(eFile);
  }
  if (arithmetic == "safe_rational") {
    using T = Rational<SafeInt64>;
    return compute_determinant_kernel<T>(eFile);
  }
  if (arithmetic == "rational") {
    using T = mpq_class;
    return compute_determinant_kernel<T>(eFile);
  }
  if (arithmetic == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return compute_determinant_kernel<T>(eFile);
  }
  if (arithmetic == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return compute_determinant_kernel<T>(eFile);
  }
  std::optional<std::string> opt_realalgebraic =
    get_postfix(arithmetic, "RealAlgebraic=");
  if (opt_realalgebraic) {
    std::string const &FileAlgebraicField = *opt_realalgebraic;
    if (!IsExistingFile(FileAlgebraicField)) {
      std::cerr << "FileAlgebraicField=" << FileAlgebraicField
                << " is missing\n";
      throw TerminalException{1};
    }
    using T_rat = mpq_class;
    HelperClassRealField<T_rat> hcrf(FileAlgebraicField);
    int const idx_real_algebraic_field = 1;
    insert_helper_real_algebraic_field(idx_real_algebraic_field, hcrf);
    using T = RealField<idx_real_algebraic_field>;
    return compute_determinant_kernel<T>(eFile);
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
    std::string eFile = argv[2];
    compute_determinant(arithmetic, eFile);
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
