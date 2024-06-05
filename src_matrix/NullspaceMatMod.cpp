// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "MAT_MatrixMod.h"
// clang-format on

template <typename T>
void full_process_type(std::string const &matrix_file, std::string const& prime,
                       std::string OutFormat, std::ostream & os_out) {
  MyMatrix<T> TheMat = ReadMatrixFile<T>(matrix_file);
  std::cerr << "TheMat=\n";
  WriteMatrix(std::cerr, TheMat);
  std::cerr << "\n";
  //
  T p = ParseScalar<T>(prime);
  //
  ResultNullspaceMod<T> res = NullspaceMatMod<T>(TheMat, p);
  //
  if (OutFormat == "GAP") {
    os_out << "return rec(dim:=" << res.dimNSP << ", NSP:=";
    WriteMatrixGAP(os_out, res.BasisTot);
    os_out << ");\n";
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat\n";
  throw TerminalException{1};
}

void process(std::string const &arith,
             std::string const &matrix_file, std::string const& prime,
             std::string OutFormat, std::ostream & os_out) {
  if (arith == "safe_rational") {
    using T = Rational<SafeInt64>;
    return full_process_type<T>(matrix_file, prime, OutFormat, os_out);
  }
  if (arith == "rational") {
    using T = mpq_class;
    return full_process_type<T>(matrix_file, prime, OutFormat, os_out);
  }
  std::cerr << "Failed to find a matching entry for arith\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 4 && argc != 6) {
      std::cerr << "This program is used as\n";
      std::cerr << "NullspaceComputation [arith] [matrix_file] [prime] [OutFormat] [OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "NullspaceComputation [arith] [matrix_file] [prime]\n";
      std::cerr << "\n";
      std::cerr << "    where\n";
      std::cerr << "arith: The arithmetic type\n";
      std::cerr << "inputMat: The input matrix\n";
      std::cerr << "prime: The string describing the prime\n";
      std::cerr << "output: The output matrix\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string matrix_file = argv[2];
    std::string prime = argv[3];
    std::string OutFormat = "GAP";
    std::string OutFile = "stderr";
    if (argc == 6) {
      OutFormat = argv[4];
      OutFile = argv[5];
    }
    //
    if (OutFile == "stderr") {
      process(arith, matrix_file, prime, OutFormat, std::cerr);
    } else {
      if (OutFile == "stdout") {
        process(arith, matrix_file, prime, OutFormat, std::cout);
      } else {
        std::ofstream os_out(OutFile);
        process(arith, matrix_file, prime, OutFormat, os_out);
      }
    }
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened in the computation\n";
    exit(e.eVal);
  }
  runtime(time1);
}
