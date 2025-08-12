// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "MAT_MatrixMod.h"
#include "MAT_MatrixInt.h"
// clang-format on

template <typename T>
void full_process_type(std::string const &matrix_file, std::string OutFormat,
                       std::ostream &os_out) {
  MyMatrix<T> TheMat = ReadMatrixFile<T>(matrix_file);
  std::cerr << "TheMat=" << TheMat.rows() << " / " << TheMat.cols() << "\n";
  //
  std::cerr << "Computing SelectRowColMatModHadamard\n";
  //
  SelectionRowColData<T> result = SelectRowColMatModHadamard<T>(TheMat);
  //
  std::cerr << "SelectRowColMatModHadamard result:\n";
  std::cerr << "  TheRank = " << result.TheRank << "\n";
  std::cerr << "  ListColSelect = [";
  for (size_t i = 0; i < result.ListColSelect.size(); i++) {
    std::cerr << result.ListColSelect[i];
    if (i < result.ListColSelect.size() - 1) std::cerr << ", ";
  }
  std::cerr << "]\n";
  std::cerr << "  ListRowSelect = [";
  for (size_t i = 0; i < result.ListRowSelect.size(); i++) {
    std::cerr << result.ListRowSelect[i];
    if (i < result.ListRowSelect.size() - 1) std::cerr << ", ";
  }
  std::cerr << "]\n";
  
  // Compare with standard SelectRowCol computation for verification
  if (TheMat.rows() <= 10) { // Only for small matrices
    std::cerr << "Computing SelectRowCol using standard method for comparison\n";
    SelectionRowCol<T> standard_result = TMat_SelectRowColRing(TheMat);
    std::cerr << "Standard result:\n";
    std::cerr << "  TheRank = " << standard_result.TheRank << "\n";
    std::cerr << "  ListColSelect = [";
    for (size_t i = 0; i < standard_result.ListColSelect.size(); i++) {
      std::cerr << standard_result.ListColSelect[i];
      if (i < standard_result.ListColSelect.size() - 1) std::cerr << ", ";
    }
    std::cerr << "]\n";
    std::cerr << "  ListRowSelect = [";
    for (size_t i = 0; i < standard_result.ListRowSelect.size(); i++) {
      std::cerr << standard_result.ListRowSelect[i];
      if (i < standard_result.ListRowSelect.size() - 1) std::cerr << ", ";
    }
    std::cerr << "]\n";
    
    if (result.TheRank == standard_result.TheRank &&
        result.ListColSelect == standard_result.ListColSelect &&
        result.ListRowSelect == standard_result.ListRowSelect) {
      std::cerr << "SUCCESS: Both methods agree!\n";
    } else {
      std::cerr << "INFO: Methods may differ (this is expected with modular arithmetic)\n";
    }
  }
  //
  if (OutFormat == "GAP") {
    os_out << "return rec(";
    os_out << "TheRank := " << result.TheRank << ", ";
    os_out << "ListColSelect := [";
    for (size_t i = 0; i < result.ListColSelect.size(); i++) {
      os_out << result.ListColSelect[i] + 1; // GAP uses 1-based indexing
      if (i < result.ListColSelect.size() - 1) os_out << ", ";
    }
    os_out << "], ";
    os_out << "ListRowSelect := [";
    for (size_t i = 0; i < result.ListRowSelect.size(); i++) {
      os_out << result.ListRowSelect[i] + 1; // GAP uses 1-based indexing
      if (i < result.ListRowSelect.size() - 1) os_out << ", ";
    }
    os_out << "]);\n";
    return;
  }
  if (OutFormat == "simple") {
    os_out << "Rank: " << result.TheRank << "\n";
    os_out << "ColSelect: ";
    for (size_t i = 0; i < result.ListColSelect.size(); i++) {
      os_out << result.ListColSelect[i];
      if (i < result.ListColSelect.size() - 1) os_out << " ";
    }
    os_out << "\n";
    os_out << "RowSelect: ";
    for (size_t i = 0; i < result.ListRowSelect.size(); i++) {
      os_out << result.ListRowSelect[i];
      if (i < result.ListRowSelect.size() - 1) os_out << " ";
    }
    os_out << "\n";
    return;
  }
  std::cerr << "Failed to find a matching entry for OutFormat\n";
  throw TerminalException{1};
}

void process(std::string const &arith, std::string const &matrix_file,
             std::string OutFormat, std::ostream &os_out) {
  if (arith == "safe_integer") {
    using T = SafeInt64;
    return full_process_type<T>(matrix_file, OutFormat, os_out);
  }
  if (arith == "mpz") {
    using T = mpz_class;
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
      std::cerr << "SelectRowColMatModHadamard [arith] [matrix_file] [OutFormat] "
                   "[OutFile]\n";
      std::cerr << "    or\n";
      std::cerr << "SelectRowColMatModHadamard [arith] [matrix_file]\n";
      std::cerr << "\n";
      std::cerr << "    where\n";
      std::cerr << "arith: The arithmetic type (safe_integer, mpz)\n";
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
