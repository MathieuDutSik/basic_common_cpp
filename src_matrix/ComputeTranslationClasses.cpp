// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheorySafeInt.h"
#include "MAT_MatrixInt.h"
// clang-format on

template <typename T>
void compute_translation_classes(std::string const &matrix_file,
                                 std::string const &OutFormat,
                                 std::ostream &os_out) {
  MyMatrix<T> M = ReadMatrixFile<T>(matrix_file);
  std::cerr << "M=" << M.rows() << " / " << M.cols() << "\n";
  std::vector<MyVector<int>> ListClasses =
      ComputeTranslationClasses<T, int>(M);
  std::cerr << "|ListClasses|=" << ListClasses.size() << "\n";
  if (OutFormat == "GAP") {
    os_out << "return [";
    for (size_t i = 0; i < ListClasses.size(); i++) {
      if (i > 0)
        os_out << ",\n";
      WriteVectorGAP(os_out, ListClasses[i]);
    }
    os_out << "];\n";
    return;
  }
  if (OutFormat == "simple") {
    for (auto &eV : ListClasses) {
      WriteVector(os_out, eV);
    }
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat\n";
  throw TerminalException{1};
}

void process(std::string const &arith, std::string const &matrix_file,
             std::string const &OutFormat, std::ostream &os_out) {
  if (arith == "safe_integer") {
    using T = SafeInt64;
    return compute_translation_classes<T>(matrix_file, OutFormat, os_out);
  }
  if (arith == "safe_rational") {
    using T = Rational<SafeInt64>;
    return compute_translation_classes<T>(matrix_file, OutFormat, os_out);
  }
  if (arith == "cpp_rational") {
    using T = boost::multiprecision::cpp_rational;
    return compute_translation_classes<T>(matrix_file, OutFormat, os_out);
  }
#ifndef DISABLE_GMP_ARITHMETIC
  if (arith == "rational") {
    using T = mpq_class;
    return compute_translation_classes<T>(matrix_file, OutFormat, os_out);
  }
  if (arith == "integer") {
    using T = mpz_class;
    return compute_translation_classes<T>(matrix_file, OutFormat, os_out);
  }
#endif
  std::cerr << "Failed to find a matching entry for arith\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "This program is used as\n";
      std::cerr << "ComputeTranslationClasses [arith] [matrix_file] "
                   "[OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "ComputeTranslationClasses [arith] [matrix_file]\n";
      std::cerr << "\n";
      std::cerr << "    where\n";
      std::cerr << "arith: The arithmetic type (safe_integer, safe_rational, "
                   "cpp_rational, rational, integer)\n";
      std::cerr << "matrix_file: The input matrix file\n";
      std::cerr << "OutFormat: GAP or simple (optional, default: simple)\n";
      std::cerr << "OutFile: Output file or stderr/stdout (optional, "
                   "default: stderr)\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string matrix_file = argv[2];
    std::string OutFormat = "simple";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    if (OutFile == "stderr") {
      process(arith, matrix_file, OutFormat, std::cerr);
    } else {
      if (OutFile == "stdout") {
        process(arith, matrix_file, OutFormat, std::cout);
      } else {
        std::ofstream os_out(OutFile);
        process(arith, matrix_file, OutFormat, os_out);
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened in the computation\n";
    exit(e.eVal);
  }
  runtime(time);
}
