// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheory.h"
#include "NumberTheorySafeInt.h"
#include "MAT_MatrixInt.h"
// clang-format on

// Tests of the fraction-free adjugate machinery:
// --- AdjugateDeterminant: A adj(A) = adj(A) A = det(A) I and the
//     determinant matching DeterminantMat, including the row-swap path
//     and the rejection of singular input.
// --- InverseFractionFreeLU on unimodular matrices.
// --- SelectIndependentRowsRing matching the greedy selection of
//     TMat_SelectRowCol over the overlying field.
// --- The ring CanonicalizeOrderedMatrix matching the computation over
//     the overlying field.

template <typename T> MyMatrix<T> RandomSquareMatrix(int n, int amp) {
  while (true) {
    MyMatrix<T> A(n, n);
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        A(i, j) = T((random() % (2 * amp + 1)) - amp);
    if (DeterminantMat(A) != 0)
      return A;
  }
}

template <typename T> void process_adjugate(int n, int nb) {
  MyMatrix<T> eId = IdentityMat<T>(n);
  for (int i = 0; i < nb; i++) {
    MyMatrix<T> A = RandomSquareMatrix<T>(n, 6);
    if (i % 3 == 1) {
      // Force the row-swap path of the elimination.
      A(0, 0) = T(0);
      if (DeterminantMat(A) == 0)
        continue;
    }
    std::pair<MyMatrix<T>, T> pair = AdjugateDeterminant(A);
    MyMatrix<T> const &adj = pair.first;
    T const &det = pair.second;
    if (det != DeterminantMat(A)) {
      std::cerr << "The determinant of AdjugateDeterminant is incorrect\n";
      std::cerr << "A=\n";
      WriteMatrix(std::cerr, A);
      throw TerminalException{1};
    }
    MyMatrix<T> prod1 = A * adj;
    MyMatrix<T> prod2 = adj * A;
    MyMatrix<T> detId = det * eId;
    if (!TestEqualityMatrix(prod1, detId) ||
        !TestEqualityMatrix(prod2, detId)) {
      std::cerr << "The relation A adj(A) = adj(A) A = det(A) I fails\n";
      std::cerr << "A=\n";
      WriteMatrix(std::cerr, A);
      throw TerminalException{1};
    }
  }
  // A singular input has to be rejected.
  MyMatrix<T> Asing = ZeroMatrix<T>(3, 3);
  Asing(0, 0) = T(1);
  Asing(1, 1) = T(1);
  bool caught = false;
  try {
    (void)AdjugateDeterminant(Asing);
  } catch (TerminalException const &e) {
    caught = true;
  }
  if (!caught) {
    std::cerr << "The singular input was not rejected\n";
    throw TerminalException{1};
  }
}

template <typename T> void process_inverse_unimodular(int n, int nb) {
  MyMatrix<T> eId = IdentityMat<T>(n);
  for (int i = 0; i < nb; i++) {
    MyMatrix<T> A = RandomUnimodularMatrix<T>(n);
    MyMatrix<T> eInv = InverseFractionFreeLU(A);
    MyMatrix<T> eProd = A * eInv;
    if (!TestEqualityMatrix(eProd, eId)) {
      std::cerr << "Error in InverseFractionFreeLU on a unimodular matrix\n";
      std::cerr << "A=\n";
      WriteMatrix(std::cerr, A);
      throw TerminalException{1};
    }
  }
}

template <typename T> MyMatrix<T> RandomRectangularMatrix(int n_row, int n_col) {
  while (true) {
    MyMatrix<T> M(n_row, n_col);
    for (int i = 0; i < n_row; i++)
      for (int j = 0; j < n_col; j++)
        M(i, j) = T((random() % 13) - 6);
    // A duplicated row and a scaled row to exercise the selection.
    if (n_row > n_col + 1) {
      M.row(n_col) = M.row(0);
      M.row(n_col + 1) = 3 * M.row(1);
    }
    if (RankMat(M) == n_col)
      return M;
  }
}

template <typename T> void process_selection(int n_row, int n_col, int nb) {
  using Tfield = typename overlying_field<T>::field_type;
  for (int i = 0; i < nb; i++) {
    MyMatrix<T> M = RandomRectangularMatrix<T>(n_row, n_col);
    std::vector<int> sel_ring = SelectIndependentRowsRing(M);
    MyMatrix<Tfield> Mf = UniversalMatrixConversion<Tfield, T>(M);
    SelectionRowCol<Tfield> eSelect = TMat_SelectRowCol(Mf);
    if (sel_ring != eSelect.ListRowSelect) {
      std::cerr << "The ring row selection differs from the field one\n";
      std::cerr << "M=\n";
      WriteMatrix(std::cerr, M);
      throw TerminalException{1};
    }
  }
}

template <typename T> void process_canonicalize(int n_row, int n_col, int nb) {
  using Tfield = typename overlying_field<T>::field_type;
  for (int i = 0; i < nb; i++) {
    MyMatrix<T> M = RandomRectangularMatrix<T>(n_row, n_col);
    MyMatrix<T> result = CanonicalizeOrderedMatrix(M);
    // The reference: the computation over the overlying field.
    MyMatrix<Tfield> Mf = UniversalMatrixConversion<Tfield, T>(M);
    MyMatrix<Tfield> OutputF = CanonicalizeOrderedMatrix_Kernel(Mf);
    MyMatrix<T> reference = UniversalMatrixConversion<T, Tfield>(OutputF);
    if (!TestEqualityMatrix(result, reference)) {
      std::cerr << "The ring CanonicalizeOrderedMatrix differs from the "
                   "field computation\n";
      std::cerr << "M=\n";
      WriteMatrix(std::cerr, M);
      std::cerr << "result=\n";
      WriteMatrix(std::cerr, result);
      std::cerr << "reference=\n";
      WriteMatrix(std::cerr, reference);
      throw TerminalException{1};
    }
  }
}

template <typename T> void process(int n) {
  int nb = 50;
  process_adjugate<T>(n, nb);
  std::cerr << "process_adjugate done\n";
  process_inverse_unimodular<T>(n, nb);
  std::cerr << "process_inverse_unimodular done\n";
  if constexpr (!is_ring_field<T>::value) {
    process_selection<T>(2 * n, n, nb);
    std::cerr << "process_selection done\n";
    process_canonicalize<T>(2 * n, n, nb);
    std::cerr << "process_canonicalize done\n";
  }
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "Test_AdjugateDeterminant [arith] [n]\n";
      return -1;
    }
    std::string arith = argv[1];
    int n;
    sscanf(argv[2], "%d", &n);
    auto f = [&]() -> void {
      if (arith == "mpz_class")
        return process<mpz_class>(n);
      if (arith == "mpq_class")
        return process<mpq_class>(n);
      if (arith == "safe_integer")
        return process<SafeInt64>(n);
      if (arith == "safe_rational")
        return process<Rational<SafeInt64>>(n);
      if (arith == "boost_cpp_int")
        return process<boost::multiprecision::cpp_int>(n);
      std::cerr << "Failed to find a matching entry\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of Test_AdjugateDeterminant\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Test_AdjugateDeterminant\n";
    exit(e.eVal);
  }
  runtime(time);
}
