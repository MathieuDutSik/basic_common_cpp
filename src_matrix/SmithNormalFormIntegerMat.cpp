// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "MAT_MatrixMod.h"
// clang-format on

template <typename T>
void full_process_type(std::string const &matrix_file, std::string OutFormat,
                       std::ostream &os_out) {
  MyMatrix<T> TheMat = ReadMatrixFile<T>(matrix_file);
  std::cerr << "TheMat=" << TheMat.rows() << " / " << TheMat.cols() << "\n";
  //
  std::cerr << "Computing Smith Normal Form\n";
  //
  MyVector<T> SNF = SmithNormalFormIntegerMat<T>(TheMat);
  //
  if (OutFormat == "GAP") {
    os_out << "return ";
    WriteVectorGAP(os_out, SNF);
    os_out << ";\n";
    return;
  }
  if (OutFormat == "simple") {
    WriteVector(os_out, SNF);
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat\n";
  throw TerminalException{1};
}

void process(std::string const &arith, std::string const &matrix_file,
             std::string OutFormat, std::ostream &os_out) {
  if (arith == "safe_rational") {
    using T = Rational<SafeInt64>;
    return full_process_type<T>(matrix_file, OutFormat, os_out);
  }
  if (arith == "mpq_class") {
    using T = mpq_class;
    return full_process_type<T>(matrix_file, OutFormat, os_out);
  }
  if (arith == "mpz_class") {
    using T = mpz_class;
    return full_process_type<T>(matrix_file, OutFormat, os_out);
  }
  if (arith == "safe_integer") {
    using T = SafeInt64;
    return full_process_type<T>(matrix_file, OutFormat, os_out);
  }
  std::cerr << "Failed to find a matching entry for arith\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3 && argc != 5) {
      std::cerr << "This program is used as\n";
      std::cerr
          << "SmithNormalFormIntegerMat [arith] [matrix_file] [OutFormat] "
             "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "SmithNormalFormIntegerMat [arith] [matrix_file]\n";
      std::cerr << "\n";
      std::cerr << "    where\n";
      std::cerr << "arith: The arithmetic type (safe_rational, mpq_class, "
                   "mpz_class, safe_integer)\n";
      std::cerr << "matrix_file: The input matrix file\n";
      std::cerr << "OutFormat: GAP or simple (optional)\n";
      std::cerr << "OutFile: Output file or stderr/stdout (optional)\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string matrix_file = argv[2];
    std::string OutFormat = "simple";
    std::string OutFile = "stderr";
    if (argc == 5) {
      OutFormat = argv[3];
      OutFile = argv[4];
    }
    //
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
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened in the computation\n";
    exit(e.eVal);
  }
  runtime(time);
}
