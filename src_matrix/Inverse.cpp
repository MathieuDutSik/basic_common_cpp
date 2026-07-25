// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheorySafeInt.h"
#include "MAT_MatrixInt.h"
// clang-format on

template <typename T>
void compute_inverse_kernel(std::string const &eFile) {
  MyMatrix<T> A = ReadMatrixFile<T>(eFile);
  int n = A.rows();
  if (n != A.cols()) {
    std::cerr << "The matrix should be square\n";
    throw TerminalException{1};
  }
  std::cerr << "Benchmarking inverse of a " << n << " x " << A.cols()
            << " matrix\n";

  MicrosecondTime time;
  MyMatrix<T> inv_dispatch = Inverse(A);
  std::cerr << "|Inverse|=" << time << "\n";

  MyMatrix<T> inv_bareiss = InverseBareiss(A);
  std::cerr << "|InverseBareiss|=" << time << "\n";

  MyMatrix<T> inv_fflu = InverseFractionFreeLU(A);
  std::cerr << "|InverseFractionFreeLU|=" << time << "\n";

  MyMatrix<T> Id = IdentityMat<T>(n);
  std::cerr << "Inverse             * A == I : " << (A * inv_dispatch == Id)
            << "\n";
  std::cerr << "InverseBareiss      * A == I : " << (A * inv_bareiss == Id)
            << "\n";
  std::cerr << "InverseFractionFreeLU * A==I : " << (A * inv_fflu == Id) << "\n";
  std::cerr << "all three agree              : "
            << (inv_dispatch == inv_bareiss && inv_bareiss == inv_fflu) << "\n";
}

void compute_inverse(std::string const &arithmetic, std::string const &eFile) {
  if (arithmetic == "safe_integer") {
    using T = SafeInt64;
    return compute_inverse_kernel<T>(eFile);
  }
  if (arithmetic == "safe_rational") {
    using T = Rational<SafeInt64>;
    return compute_inverse_kernel<T>(eFile);
  }
  if (arithmetic == "cpp_rational") {
    using T = boost::multiprecision::cpp_rational;
    return compute_inverse_kernel<T>(eFile);
  }
#ifndef DISABLE_GMP_ARITHMETIC
  if (arithmetic == "integer") {
    using T = mpz_class;
    return compute_inverse_kernel<T>(eFile);
  }
  if (arithmetic == "rational") {
    using T = mpq_class;
    return compute_inverse_kernel<T>(eFile);
  }
  if (arithmetic == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return compute_inverse_kernel<T>(eFile);
  }
  if (arithmetic == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return compute_inverse_kernel<T>(eFile);
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
    return compute_inverse_kernel<T>(eFile);
  }
#endif
  std::cerr << "Failed to find a matching arithmetic\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "Inverse [arithmetic] [inputMat]\n";
      return -1;
    }
    std::string arithmetic = argv[1];
    std::string eFile = argv[2];
    compute_inverse(arithmetic, eFile);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
