// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "QuadField.h"
#include "MAT_Matrix.h"
// clang-format on

template<typename T>
void full_process_type(std::string const& input, std::string const& output) {
  MyMatrix<T> TheMat = ReadMatrixFile<T>(input);
  std::cerr << "TheMat=\n";
  WriteMatrix(std::cerr, TheMat);
  std::cerr << "\n";
  //
  MyMatrix<T> TheKer = NullspaceMat(TheMat);
  std::cerr << "TheKer=\n";
  WriteMatrix(std::cerr, TheKer);
  std::cerr << "\n";
  //
  std::ofstream os(output);
  os << "return ";
  WriteMatrixGAP(os, TheKer);
  os << ";\n";
}


void process(std::string const& arith, std::string const& input, std::string const& output) {
  if (arith == "safe_integer") {
    using T = SafeInt64;
    return full_process_type<T>(input, output);
  }
  if (arith == "safe_rational") {
    using T = Rational<SafeInt64>;
    return full_process_type<T>(input, output);
  }
  if (arith == "rational") {
    using T = mpq_class;
    return full_process_type<T>(input, output);
  }
  if (arith == "Qsqrt5") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 5>;
    return full_process_type<T>(input, output);
  }
  if (arith == "Qsqrt2") {
    using Trat = mpq_class;
    using T = QuadField<Trat, 2>;
    return full_process_type<T>(input, output);
  }
  std::optional<std::string> opt_realalgebraic =
    get_postfix(arith, "RealAlgebraic=");
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
    return full_process_type<T>(input, output);
  }
  std::cerr << "Failed to find a matching entry for arith\n";
  throw TerminalException{1};
}



int main(int argc, char *argv[]) {
  SingletonTime time1;
  try {
    if (argc != 4) {
      fprintf(stderr, "Number of argument is = %d\n", argc);
      fprintf(stderr, "This program is used as\n");
      fprintf(stderr, "NullspaceComputation [arith] [inputMat] [output]\n");
      fprintf(stderr, "    where\n");
      fprintf(stderr, "arith: The arithmetic type\n");
      fprintf(stderr, "inputMat: The input matrix\n");
      fprintf(stderr, "output: The output matrix\n");
      return -1;
    }
    //
    std::string arith = argv[1];
    std::string input = argv[2];
    std::string output = argv[3];
    //
    process(arith, input, output);
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened in the computation\n";
    exit(e.eVal);
  }
  runtime(time1);
}
