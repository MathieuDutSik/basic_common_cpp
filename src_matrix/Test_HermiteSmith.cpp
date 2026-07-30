// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheory.h"
#include "NumberTheorySafeInt.h"
#include "MAT_MatrixInt.h"
// clang-format on

// Tests of the Hermite and Smith normal forms:
// --- ComputeRowHermiteNormalForm: H = U M with U unimodular, echelon
//     shape, and canonicity (invariance under left multiplication by a
//     unimodular matrix, idempotence).
// --- ComputeColHermiteNormalForm: the transposed statements.
// --- SmithNormalForm: ROW M COL diagonal with unimodular ROW and COL,
//     the divisibility chain of the invariant factors, consistency with
//     SmithNormalFormInvariant, and invariance of the invariant factors
//     under unimodular multiplication on both sides.

template <typename T> MyMatrix<T> RandomMatrix(int n_row, int n_col, int amp) {
  MyMatrix<T> A(n_row, n_col);
  for (int i = 0; i < n_row; i++)
    for (int j = 0; j < n_col; j++)
      A(i, j) = T((random() % (2 * amp + 1)) - amp);
  return A;
}

// A mild unimodular matrix (few shears with small pivots), so that the
// conjugated matrices stay within the range of SafeInt64.
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

template <typename T> bool IsUnimodular(MyMatrix<T> const &U) {
  T det = DeterminantMat(U);
  return det == 1 || det == -1;
}

// The row echelon shape: the leading column indices are strictly
// increasing and the zero rows are at the bottom.
template <typename T> bool IsRowEchelon(MyMatrix<T> const &H) {
  int n_row = H.rows();
  int n_col = H.cols();
  int prev = -1;
  bool zero_seen = false;
  for (int i = 0; i < n_row; i++) {
    int lead = -1;
    for (int j = 0; j < n_col; j++)
      if (H(i, j) != 0) {
        lead = j;
        break;
      }
    if (lead == -1) {
      zero_seen = true;
    } else {
      if (zero_seen || lead <= prev)
        return false;
      prev = lead;
    }
  }
  return true;
}

template <typename T> void process_hnf(int n, int nb) {
  for (int i = 0; i < nb; i++) {
    int n_row = n + (i % 3);
    MyMatrix<T> M = RandomMatrix<T>(n_row, n, 6);
    std::pair<MyMatrix<T>, MyMatrix<T>> pairRow = ComputeRowHermiteNormalForm(M);
    MyMatrix<T> const &U = pairRow.first;
    MyMatrix<T> const &H = pairRow.second;
    MyMatrix<T> prod = U * M;
    if (!TestEqualityMatrix(prod, H) || !IsUnimodular(U) || !IsRowEchelon(H)) {
      std::cerr << "The row Hermite normal form H = U M fails\n";
      std::cerr << "M=\n";
      WriteMatrix(std::cerr, M);
      throw TerminalException{1};
    }
    // Canonicity: a unimodular change of generators does not change H.
    MyMatrix<T> V = MildUnimodularMatrix<T>(n_row);
    MyMatrix<T> M2 = V * M;
    MyMatrix<T> H2 = ComputeRowHermiteNormalForm(M2).second;
    MyMatrix<T> H3 = ComputeRowHermiteNormalForm(H).second;
    if (!TestEqualityMatrix(H2, H) || !TestEqualityMatrix(H3, H)) {
      std::cerr << "The row Hermite normal form is not canonical\n";
      std::cerr << "M=\n";
      WriteMatrix(std::cerr, M);
      throw TerminalException{1};
    }
    // The column version on the transposed statements.
    MyMatrix<T> Mtr = TransposedMat(M);
    std::pair<MyMatrix<T>, MyMatrix<T>> pairCol =
        ComputeColHermiteNormalForm(Mtr);
    MyMatrix<T> prodCol = Mtr * pairCol.first;
    if (!TestEqualityMatrix(prodCol, pairCol.second) ||
        !IsUnimodular(pairCol.first)) {
      std::cerr << "The column Hermite normal form H = M U fails\n";
      std::cerr << "Mtr=\n";
      WriteMatrix(std::cerr, Mtr);
      throw TerminalException{1};
    }
    MyMatrix<T> Hcol_tr = TransposedMat(pairCol.second);
    if (!IsRowEchelon(Hcol_tr)) {
      std::cerr << "The column Hermite normal form is not column echelon\n";
      std::cerr << "Mtr=\n";
      WriteMatrix(std::cerr, Mtr);
      throw TerminalException{1};
    }
  }
}

template <typename T> void process_smith(int n, int nb) {
  for (int i = 0; i < nb; i++) {
    int n_row = n + (i % 3);
    MyMatrix<T> M = RandomMatrix<T>(n_row, n, 6);
    if (i % 3 == 2) {
      // A rank deficient input: last row a combination of the others.
      MyVector<T> comb = ZeroVector<T>(n);
      for (int k = 0; k + 1 < n_row; k++) {
        T c = T((random() % 5) - 2);
        for (int j = 0; j < n; j++)
          comb(j) += c * M(k, j);
      }
      for (int j = 0; j < n; j++)
        M(n_row - 1, j) = comb(j);
    }
    ResultSmithNormalForm<T> res = SmithNormalForm(M);
    if (!IsUnimodular(res.ROW) || !IsUnimodular(res.COL)) {
      std::cerr << "The Smith normal form transformations are not unimodular\n";
      throw TerminalException{1};
    }
    MyMatrix<T> P = res.ROW * M * res.COL;
    int len = res.Invariant.size();
    for (int i_row = 0; i_row < P.rows(); i_row++) {
      for (int i_col = 0; i_col < P.cols(); i_col++) {
        T expect(0);
        if (i_row == i_col && i_row < len)
          expect = res.Invariant(i_row);
        if (P(i_row, i_col) != expect) {
          std::cerr << "ROW M COL does not match the invariant factors\n";
          std::cerr << "M=\n";
          WriteMatrix(std::cerr, M);
          std::cerr << "P=\n";
          WriteMatrix(std::cerr, P);
          throw TerminalException{1};
        }
      }
    }
    // The divisibility chain, with the zero factors at the end.
    for (int k = 0; k + 1 < len; k++) {
      T const &a = res.Invariant(k);
      T const &b = res.Invariant(k + 1);
      bool ok = [&]() -> bool {
        if (a == 0)
          return b == 0;
        return ResInt(b, a) == 0;
      }();
      if (!ok) {
        std::cerr << "The divisibility chain of the invariant factors fails\n";
        std::cerr << "M=\n";
        WriteMatrix(std::cerr, M);
        throw TerminalException{1};
      }
    }
    // Consistency with SmithNormalFormInvariant and invariance under
    // unimodular multiplication on both sides.
    MyVector<T> inv_direct = SmithNormalFormInvariant(M);
    MyMatrix<T> M2 =
        MildUnimodularMatrix<T>(n_row) * M * MildUnimodularMatrix<T>(n);
    MyVector<T> inv_conj = SmithNormalFormInvariant(M2);
    auto abs_equal = [&](MyVector<T> const &V1, MyVector<T> const &V2) -> bool {
      if (V1.size() != V2.size())
        return false;
      for (int k = 0; k < V1.size(); k++)
        if (T_abs(V1(k)) != T_abs(V2(k)))
          return false;
      return true;
    };
    if (!abs_equal(inv_direct, res.Invariant) ||
        !abs_equal(inv_conj, res.Invariant)) {
      std::cerr << "The Smith invariant factors are not invariant\n";
      std::cerr << "M=\n";
      WriteMatrix(std::cerr, M);
      std::cerr << "res.Invariant=";
      WriteVector(std::cerr, res.Invariant);
      std::cerr << "inv_direct=";
      WriteVector(std::cerr, inv_direct);
      std::cerr << "inv_conj=";
      WriteVector(std::cerr, inv_conj);
      throw TerminalException{1};
    }
  }
}

template <typename T> void process(int n) {
  int nb = 30;
  process_hnf<T>(n, nb);
  std::cerr << "process_hnf done\n";
  process_smith<T>(n, nb);
  std::cerr << "process_smith done\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "Test_HermiteSmith [arith] [n]\n";
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
      std::cerr << "Failed to find a matching entry\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of Test_HermiteSmith\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Test_HermiteSmith\n";
    exit(e.eVal);
  }
  runtime(time);
}
