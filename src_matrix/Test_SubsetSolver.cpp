// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheorySafeInt.h"
#include "Boost_bitset.h"
#include "MAT_Matrix_SubsetSolver.h"
#include "MAT_MatrixRankmat.h"
// clang-format on

// Tests of SubsetRankOneSolver: for a subset of d-1 independent rows of
// a m x d matrix, the kernel vector is non-zero and orthogonal to the
// selected rows, and the accelerated implementation agrees with the
// field one up to a positive scalar (GetPositiveKernelVector fixes the
// sign against the first row outside the subset).

template <typename T> MyMatrix<T> RandomIntegralMatrix(int n_row, int n_col) {
  MyMatrix<T> A(n_row, n_col);
  for (int i = 0; i < n_row; i++)
    for (int j = 0; j < n_col; j++)
      A(i, j) = T((random() % 13) - 6);
  return A;
}

// A random subset of n_col - 1 rows spanning a rank n_col - 1 subspace,
// whose first complement row is not in that span (the sign convention of
// GetPositiveKernelVector needs a non-zero scalar product there).
template <typename T> Face RandomCoRankOneFace(MyMatrix<T> const &EXT) {
  int n_row = EXT.rows();
  int n_col = EXT.cols();
  while (true) {
    std::vector<int> idx(n_row);
    for (int i = 0; i < n_row; i++)
      idx[i] = i;
    for (int i = n_row - 1; i > 0; i--) {
      int j = random() % (i + 1);
      std::swap(idx[i], idx[j]);
    }
    Face f(n_row);
    for (int k = 0; k < n_col - 1; k++)
      f[idx[k]] = 1;
    std::vector<int> selected;
    for (int i = 0; i < n_row; i++)
      if (f[i])
        selected.push_back(i);
    MyMatrix<T> Msel = SelectRow(EXT, selected);
    if (RankMat(Msel) != n_col - 1)
      continue;
    int i_first = -1;
    for (int i = 0; i < n_row; i++)
      if (f[i] == 0) {
        i_first = i;
        break;
      }
    MyMatrix<T> Mext(n_col, n_col);
    for (int k = 0; k + 1 < n_col; k++)
      Mext.row(k) = Msel.row(k);
    Mext.row(n_col - 1) = EXT.row(i_first);
    if (RankMat(Mext) != n_col)
      continue;
    return f;
  }
}

template <typename T> void process(int n) {
  using Tint = typename SubsetRankOneSolver<T>::Tint;
  int nb = 50;
  for (int i = 0; i < nb; i++) {
    int n_row = n + 3 + (i % 3);
    MyMatrix<T> EXT = RandomIntegralMatrix<T>(n_row, n);
    if (RankMat(EXT) != n)
      continue;
    MyMatrix<Tint> EXT_int = UniversalMatrixConversion<Tint, T>(EXT);
    Face f = RandomCoRankOneFace(EXT);
    SubsetRankOneSolver<T> solver(EXT_int);
    SubsetRankOneSolver_Field<T> solver_field(EXT);
    MyVector<Tint> V1 = solver.GetPositiveKernelVector(f);
    MyVector<T> V2 = solver_field.GetPositiveKernelVector(f);
    if (IsZeroVector(V1) || IsZeroVector(V2)) {
      std::cerr << "The kernel vector is zero\n";
      throw TerminalException{1};
    }
    // Orthogonality to the selected rows.
    for (int iRow = 0; iRow < n_row; iRow++) {
      if (f[iRow] == 1) {
        Tint scal(0);
        for (int iCol = 0; iCol < n; iCol++)
          scal += EXT_int(iRow, iCol) * V1(iCol);
        if (scal != 0) {
          std::cerr << "The kernel vector is not orthogonal to the subset\n";
          WriteMatrix(std::cerr, EXT);
          throw TerminalException{1};
        }
      }
    }
    // The two implementations agree up to a positive scalar.
    MyVector<T> V1_T = UniversalVectorConversion<T, Tint>(V1);
    for (int iCol = 0; iCol < n; iCol++) {
      for (int jCol = iCol + 1; jCol < n; jCol++) {
        if (V1_T(iCol) * V2(jCol) != V1_T(jCol) * V2(iCol)) {
          std::cerr << "The two kernel vectors are not proportional\n";
          WriteMatrix(std::cerr, EXT);
          throw TerminalException{1};
        }
      }
    }
    for (int iCol = 0; iCol < n; iCol++) {
      if ((V1_T(iCol) > 0 && V2(iCol) < 0) ||
          (V1_T(iCol) < 0 && V2(iCol) > 0)) {
        std::cerr << "The signs of the two kernel vectors differ\n";
        WriteMatrix(std::cerr, EXT);
        throw TerminalException{1};
      }
    }
  }
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "Test_SubsetSolver [arith] [n]\n";
      return -1;
    }
    std::string arith = argv[1];
    int n;
    sscanf(argv[2], "%d", &n);
    auto f = [&]() -> void {
      if (arith == "mpq_class")
        return process<mpq_class>(n);
      if (arith == "safe_rational")
        return process<Rational<SafeInt64>>(n);
      std::cerr << "Failed to find a matching entry\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of Test_SubsetSolver\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Test_SubsetSolver\n";
    exit(e.eVal);
  }
  runtime(time);
}
