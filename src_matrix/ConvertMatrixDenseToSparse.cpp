// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// Rewrites a matrix from the dense text format
//
//     nbRow nbCol
//     <nbRow * nbCol values>
//
// into the sparse one read by ReadSparseMatrix
//
//     nbRow nbCol nnz
//     iRow iCol value
//     ...
//
// The entries are streamed rather than held in a matrix, so the memory is
// proportional to the number of NONZERO entries and a matrix far too
// large to hold densely can still be converted. That is the case that
// motivates the program: a 11568 x 12119 boundary matrix of a chain
// complex is 280 MB of text and 2.3 GB once held densely, for 81720
// entries that actually carry a value.
//
// clang-format off
#include "NumberTheory.h"
#include "MAT_SparseMatrix.h"
#include "Timings.h"
#include <fstream>
#include <string>
#include <vector>
// clang-format on

template <typename T>
void process(std::string const &FileI, std::string const &FileO) {
  std::ifstream is(FileI);
  if (!is) {
    std::cerr << "Failed to open " << FileI << "\n";
    throw TerminalException{1};
  }
  int nbRow, nbCol;
  is >> nbRow >> nbCol;
  if (!is) {
    std::cerr << "Failed to read the dimensions from " << FileI << "\n";
    throw TerminalException{1};
  }
  struct Entry {
    int iRow;
    int iCol;
    T val;
  };
  std::vector<Entry> entries;
  T val;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    for (int iCol = 0; iCol < nbCol; iCol++) {
      is >> val;
      if (!is) {
        std::cerr << "Failed to read the entry (" << iRow << "," << iCol
                  << ") of " << FileI << "\n";
        throw TerminalException{1};
      }
      if (val != 0)
        entries.push_back({iRow, iCol, val});
    }
  }
  std::ofstream os(FileO);
  os << nbRow << " " << nbCol << " " << entries.size() << "\n";
  for (auto &ent : entries)
    os << ent.iRow << " " << ent.iCol << " " << ent.val << "\n";
  std::cerr << "|M|=" << nbRow << " / " << nbCol << " nnz=" << entries.size()
            << " density="
            << 100.0 * static_cast<double>(entries.size()) /
                   (static_cast<double>(nbRow) * nbCol)
            << "%\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 4) {
      std::cerr << "This program is used as\n";
      std::cerr << "ConvertMatrixDenseToSparse [arith] [inputMat] [output]\n";
      std::cerr << "---\n";
      std::cerr << "arith: mpz_class, mpq_class, safe_integer\n";
      std::cerr << "inputMat: the matrix in the dense format\n";
      std::cerr << "output: the matrix in the sparse format\n";
      return -1;
    }
    std::string arith = argv[1];
    std::string FileI = argv[2];
    std::string FileO = argv[3];
    auto f = [&]() -> void {
      if (arith == "mpz_class")
        return process<mpz_class>(FileI, FileO);
      if (arith == "mpq_class")
        return process<mpq_class>(FileI, FileO);
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
