// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheory.h"
#include "NumberTheorySafeInt.h"
#include "MAT_MatrixInt.h"
// clang-format on

// Tests of the rank computation:
// --- The ring RankMat (fraction-free Bareiss) against the field
//     Gaussian elimination over the overlying field.
// --- Planted low rank factorizations M = C B.
// --- Invariance under transposition, mild unimodular multiplication on
//     both sides, appended row combinations and appended zero rows.
// --- The edge cases of the zero matrix and the identity.

template <typename T> MyMatrix<T> RandomMatrix(int n_row, int n_col, int amp) {
  MyMatrix<T> A(n_row, n_col);
  for (int i = 0; i < n_row; i++)
    for (int j = 0; j < n_col; j++)
      A(i, j) = T((random() % (2 * amp + 1)) - amp);
  return A;
}

template <typename T> MyMatrix<T> MildUnimodularMatrix(int n) {
  MyMatrix<T> RetMat = IdentityMat<T>(n);
  for (int iter = 0; iter < n; iter++) {
    int idx1 = random() % n;
    int idx2 = random() % n;
    if (idx1 != idx2) {
      MyMatrix<T> eMat = IdentityMat<T>(n);
      eMat(idx1, idx2) = T((random() % 5) - 2);
      RetMat = eMat * RetMat;
    }
  }
  return RetMat;
}

// The independent reference: Gaussian elimination over the overlying
// field (the field code path of RankMat).
template <typename T> int RankReference(MyMatrix<T> const &M) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> Mf = UniversalMatrixConversion<Tfield, T>(M);
  return RankMatKernel(Mf);
}

template <typename T> void process(int n) {
  int nb = 100;
  for (int i = 0; i < nb; i++) {
    // A planted factorization M = C B of rank at most r.
    int n_row = n + (i % 3);
    int n_col = n + 1;
    int r = 1 + (random() % n);
    MyMatrix<T> C = RandomMatrix<T>(n_row, r, 4);
    MyMatrix<T> B = RandomMatrix<T>(r, n_col, 4);
    MyMatrix<T> M = C * B;
    int rank = RankMat(M);
    int rank_ref = RankReference(M);
    if (rank != rank_ref || rank > r) {
      std::cerr << "RankMat disagrees with the field reference\n";
      std::cerr << "rank=" << rank << " rank_ref=" << rank_ref << " r=" << r
                << "\n";
      std::cerr << "M=\n";
      WriteMatrix(std::cerr, M);
      throw TerminalException{1};
    }
    // Transposition.
    MyMatrix<T> Mtr = TransposedMat(M);
    if (RankMat(Mtr) != rank) {
      std::cerr << "The rank is not invariant under transposition\n";
      WriteMatrix(std::cerr, M);
      throw TerminalException{1};
    }
    // Mild unimodular multiplication on both sides.
    MyMatrix<T> M2 = MildUnimodularMatrix<T>(n_row) * M * MildUnimodularMatrix<T>(n_col);
    if (RankMat(M2) != rank) {
      std::cerr << "The rank is not invariant under unimodular equivalence\n";
      WriteMatrix(std::cerr, M);
      throw TerminalException{1};
    }
    // An appended row combination and an appended zero row do not change
    // the rank.
    MyMatrix<T> M3(n_row + 2, n_col);
    for (int k = 0; k < n_row; k++)
      M3.row(k) = M.row(k);
    M3.row(n_row) = M.row(0) + T(2) * M.row(n_row - 1);
    for (int j = 0; j < n_col; j++)
      M3(n_row + 1, j) = T(0);
    if (RankMat(M3) != rank) {
      std::cerr << "The rank changed under appended dependent rows\n";
      WriteMatrix(std::cerr, M);
      throw TerminalException{1};
    }
  }
  // The edge cases.
  MyMatrix<T> Zer = ZeroMatrix<T>(n, n + 2);
  if (RankMat(Zer) != 0) {
    std::cerr << "The zero matrix has non-zero rank\n";
    throw TerminalException{1};
  }
  MyMatrix<T> Id = IdentityMat<T>(n);
  if (RankMat(Id) != n) {
    std::cerr << "The identity matrix has the wrong rank\n";
    throw TerminalException{1};
  }
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "Test_RankMat [arith] [n]\n";
      return -1;
    }
    std::string arith = argv[1];
    int n;
    sscanf(argv[2], "%d", &n);
    auto f = [&]() -> void {
      if (arith == "mpz_class")
        return process<mpz_class>(n);
      if (arith == "safe_integer")
        return process<SafeInt64>(n);
      if (arith == "boost_cpp_int")
        return process<boost::multiprecision::cpp_int>(n);
      if (arith == "mpq_class")
        return process<mpq_class>(n);
      if (arith == "safe_rational")
        return process<Rational<SafeInt64>>(n);
      std::cerr << "Failed to find a matching entry\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of Test_RankMat\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Test_RankMat\n";
    exit(e.eVal);
  }
  runtime(time);
}
