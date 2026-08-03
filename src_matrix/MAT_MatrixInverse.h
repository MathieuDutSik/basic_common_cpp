// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXINVERSE_H_
#define SRC_MATRIX_MAT_MATRIXINVERSE_H_

// Matrix inverse over fields and rings.
//
//  --- TMat_Inverse_destroy: in-place Gauss-Jordan core (destroys its input).
//  --- InverseKernel / Inverse_destroy: thin wrappers over the core.
//  --- Inverse: dispatch entry point (field directly, other rings via the
//      overlying field).
//  --- CongrMap: the congruence map M -> (M^-1)^T.

// clang-format off
#include "MAT_MatrixFund.h"
// clang-format on

#ifdef SANITY_CHECK
#define SANITY_CHECK_MATRIX_INVERSE
#endif

template <typename T>
void TMat_Inverse_destroy(MyMatrix<T> &Input, MyMatrix<T> &Output) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in TMat_Inverse_destroy");
  int iCol, iRow;
  int iRowB;
  int nbRow = Input.rows();
  int nbCol = Input.cols();
  T prov1;
#ifdef DEBUG_MAT_MATRIX_DISABLE
  std::cerr << "TMat_Inverse_destroy, step 1\n";
#endif
#ifdef SANITY_CHECK_MAT_MATRIX
  if (nbRow != nbCol) {
    std::cerr << "Error on nbRow, nbCol in TMat_Inverse_destroy";
    throw TerminalException{1};
  }
#endif
  for (iCol = 0; iCol < nbRow; iCol++)
    for (iRow = 0; iRow < nbRow; iRow++) {
      if (iRow == iCol)
        prov1 = 1;
      else
        prov1 = 0;
      Output(iRow, iCol) = prov1;
    }
#ifdef DEBUG_MAT_MATRIX_DISABLE
  std::cerr << "TMat_Inverse_destroy, step 2\n";
#endif
  int iColFound;
  for (iRow = 0; iRow < nbRow; iRow++) {
#ifdef DEBUG_MAT_MATRIX_DISABLE
    std::cerr << "iRow=" << iRow << "\n";
    std::cerr << "Input=\n";
    WriteMatrix(std::cerr, Input);
#endif
    iColFound = SelectBestPivot<T>(
        nbCol, [&](int iCol) -> T const & { return Input(iRow, iCol); },
        [&](int iCol) -> bool { return iCol >= iRow; });
#ifdef SANITY_CHECK_MAT_MATRIX
    if (iColFound == -1) {
      std::cerr << "Error during the computation of the matrix inverse\n";
      throw TerminalException{1};
    }
#endif
    prov1 = 1 / Input(iRow, iColFound);
    for (iRowB = 0; iRowB < nbRow; iRowB++)
      Output(iRowB, iColFound) *= prov1;
    for (iRowB = iRow; iRowB < nbRow; iRowB++)
      Input(iRowB, iColFound) *= prov1;
    for (iCol = 0; iCol < nbCol; iCol++)
      if (iCol != iColFound) {
        prov1 = Input(iRow, iCol);
        if (prov1 != 0) {
          // Y -= prov1 * pivot_col. For types where the product allocates a
          // temporary each evaluation (mpz_class, mpq_class, ...) reuse prov2;
          // otherwise the fused form is at least as good (see is_fma_prefered).
          for (iRowB = 0; iRowB < nbRow; iRowB++)
            SubMul(Output(iRowB, iCol), prov1, Output(iRowB, iColFound));
          for (iRowB = iRow; iRowB < nbRow; iRowB++)
            SubMul(Input(iRowB, iCol), prov1, Input(iRowB, iColFound));
        }
      }
    if (iColFound != iRow) {
      for (iRowB = 0; iRowB < nbRow; iRowB++)
        std::swap(Output(iRowB, iColFound), Output(iRowB, iRow));
      for (iRowB = iRow; iRowB < nbRow; iRowB++)
        std::swap(Input(iRowB, iColFound), Input(iRowB, iRow));
    }
  }
#ifdef DEBUG_MAT_MATRIX_DISABLE
  std::cerr << "TMat_Inverse_destroy, step 3\n";
#endif
}

template <typename T> MyMatrix<T> InverseKernel(MyMatrix<T> const &Input) {
  int nbRow = Input.rows();
  MyMatrix<T> provMat = Input;
  MyMatrix<T> Output(nbRow, nbRow);
  TMat_Inverse_destroy(provMat, Output);
  return Output;
}

template <typename T> MyMatrix<T> Inverse_destroy(MyMatrix<T> &Input) {
  int nbRow = Input.rows();
  MyMatrix<T> Output(nbRow, nbRow);
  TMat_Inverse_destroy(Input, Output);
  return Output;
}

// Fraction-free Gauss-Jordan elimination, also known as fraction-free
// Gauss-Jordan or the Bareiss-Montante algorithm. Gauss-Jordan is run on the
// augmented matrix [A | I] using only ring operations and EXACT divisions: at
// step k every entry (pivot * M(i,j) - M(i,k) * M(k,j)) is divisible by the
// previous pivot, because Bareiss's theorem makes each intermediate entry a
// minor determinant of the input. This drives [A | I] to [det(A) * I | adj(A)],
// where adj(A) is the adjugate, so A^{-1} = adj(A) / det(A). Keeping the
// elimination fraction-free bounds intermediate operand growth (Hadamard), so
// for exact heavy arithmetic (rationals, number fields) it is faster than the
// classical Gauss-Jordan of InverseKernel, in the same way that
// DeterminantMatBareiss beats plain Gaussian elimination.
//
// Valid over any integral domain. The final division adj(A) / det(A) is exact
// over a field; over a ring of integers it is exact only when A is unimodular,
// and a non-exact division is reported rather than silently truncated.
template <typename T> MyMatrix<T> InverseBareiss(MyMatrix<T> const &Input) {
  int n = Input.rows();
  // Augmented matrix M = [A | I], of size n x 2n.
  MyMatrix<T> M(n, 2 * n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      M(i, j) = Input(i, j);
      M(i, n + j) = (i == j) ? T(1) : T(0);
    }
  T prev(1);
  for (int k = 0; k < n; k++) {
    if (M(k, k) == 0) {
      int r = -1;
      for (int i = k + 1; i < n; i++)
        if (M(i, k) != 0) {
          r = i;
          break;
        }
      if (r == -1) {
        std::cerr << "InverseBareiss: the matrix is singular\n";
        throw TerminalException{1};
      }
      M.row(k).swap(M.row(r));
    }
    T pivot = M(k, k);
    for (int i = 0; i < n; i++) {
      if (i == k)
        continue;
      for (int j = 0; j < 2 * n; j++) {
        if (j == k)
          continue;
        T val = pivot * M(i, j) - M(i, k) * M(k, j);
        T quot = val / prev; // exact division guaranteed by Bareiss's theorem
#ifdef DEBUG_MAT_MATRIX
        if (quot * prev != val) {
          std::cerr << "InverseBareiss: non-exact division, T is not an "
                       "integral domain\n";
          throw TerminalException{1};
        }
#endif
        M(i, j) = quot;
      }
      M(i, k) = T(0);
    }
    prev = pivot;
  }
  // Left block is now det(A) * I (the shared diagonal value) and the right
  // block is the adjugate adj(A) = det(A) * A^{-1}; the swap sign cancels in the
  // ratio, so no sign bookkeeping is needed.
  T det = prev;
  MyMatrix<T> Output(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      T num = M(i, n + j);
      T quot = num / det;
      if (quot * det != num) {
        std::cerr << "InverseBareiss: A^{-1} is not representable over T "
                     "(the determinant does not divide the adjugate)\n";
        throw TerminalException{1};
      }
      Output(i, j) = quot;
    }
  return Output;
}

// Fraction-free LU factorization (Zhou & Jeffrey, "Fraction-free matrix factors:
// new forms for LU and QR factors", 2008). It writes
//     A = L * D^{-1} * U
// with integer lower-triangular L, integer upper-triangular U, and diagonal
//     D = diag(d_{-1} d_0, d_0 d_1, ..., d_{n-2} d_{n-1}),  d_{-1} = 1,
// where d_k is the pivot at step k, i.e. the leading (k+1) x (k+1) minor of A;
// det(A) = d_{n-1}. L and U carry the pivots d_k on their diagonals, every entry
// is a subdeterminant of A, and all divisions are exact (Bareiss's theorem). No
// pivoting is done, so this form requires the leading principal minors to be
// non-zero (a generic dense matrix qualifies).
template <typename T> struct FractionFreeLU_result {
  MyMatrix<T> L;
  MyMatrix<T> U;
  std::vector<T> d; // pivots d_0 .. d_{n-1}; det(A) = d_{n-1}
};

template <typename T>
FractionFreeLU_result<T> FractionFreeLU(MyMatrix<T> const &A) {
  int n = A.rows();
  MyMatrix<T> M = A;
  MyMatrix<T> L = ZeroMatrix<T>(n, n);
  std::vector<T> d(n);
  T prev(1);
  for (int k = 0; k < n; k++) {
    // Column k of L is the current column-k below (and on) the diagonal.
    for (int i = k; i < n; i++)
      L(i, k) = M(i, k);
    T pivot = M(k, k);
    d[k] = pivot;
    for (int i = k + 1; i < n; i++) {
      for (int j = k + 1; j < n; j++) {
        T val = pivot * M(i, j) - M(i, k) * M(k, j);
        M(i, j) = val / prev; // exact division (Bareiss)
      }
      M(i, k) = T(0);
    }
    prev = pivot;
  }
  MyMatrix<T> U = ZeroMatrix<T>(n, n);
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      U(i, j) = M(i, j);
  return {std::move(L), std::move(U), std::move(d)};
}

// The adjugate matrix and the determinant, computed fraction-free through
// the LU approach. Rather than reducing all the way to a diagonal like the
// Bareiss-Montante Gauss-Jordan (InverseBareiss), it performs fraction-free
// forward elimination on [A | I] only below the diagonal -- yielding the
// upper-triangular U on the left and the forward-transformed identity C on
// the right -- and then a fraction-free back substitution, one column at a
// time. For column c the scaled solution
//     xhat_i = (det * C(i,c) - sum_{k>i} U(i,k) xhat_k) / d_i
// stays integral (Bareiss), and xhat is the corresponding column of the
// adjugate adj(A). The returned pair (adj(A), det(A)) satisfies
//     A * adj(A) = adj(A) * A = det(A) * I_n
// with every entry of adj(A) in the ring (a signed cofactor of A), so the
// whole computation runs over a ring (only exact divisions occur).
// Precondition: A is non-singular (the elimination needs a non-zero pivot
// in every column; singular input throws).
template <typename T>
std::pair<MyMatrix<T>, T> AdjugateDeterminant(MyMatrix<T> const &Input) {
  int n = Input.rows();
  // Augmented matrix M = [A | I], of size n x 2n.
  MyMatrix<T> M(n, 2 * n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      M(i, j) = Input(i, j);
      M(i, n + j) = (i == j) ? T(1) : T(0);
    }
  T prev(1);
  bool neg = false;
  // Fraction-free forward elimination (only the rows below the pivot).
  for (int k = 0; k < n; k++) {
    if (M(k, k) == 0) {
      int r = -1;
      for (int i = k + 1; i < n; i++)
        if (M(i, k) != 0) {
          r = i;
          break;
        }
      if (r == -1) {
        std::cerr << "AdjugateDeterminant: the matrix is singular\n";
        throw TerminalException{1};
      }
      M.row(k).swap(M.row(r));
      neg = !neg;
    }
    T pivot = M(k, k);
    for (int i = k + 1; i < n; i++) {
      for (int j = k + 1; j < 2 * n; j++) {
        T val = pivot * M(i, j) - M(i, k) * M(k, j);
        T quot = val / prev; // exact division (Bareiss)
#ifdef SANITY_CHECK_MATRIX_INVERSE
        if (quot * prev != val) {
          std::cerr << "AdjugateDeterminant: non-exact elimination division\n";
          throw TerminalException{1};
        }
#endif
        M(i, j) = quot;
      }
      M(i, k) = T(0);
    }
    prev = pivot;
  }
  T det = neg ? -prev : prev;
  // Fraction-free back substitution, one column of the adjugate at a time.
  MyMatrix<T> Output(n, n);
  std::vector<T> xhat(n);
  for (int c = 0; c < n; c++) {
    for (int i = n - 1; i >= 0; i--) {
      T s = det * M(i, n + c);
      for (int k = i + 1; k < n; k++)
        SubMul(s, M(i, k), xhat[k]);
      T quot = s / M(i, i); // exact division (Bareiss)
#ifdef SANITY_CHECK_MATRIX_INVERSE
      if (quot * M(i, i) != s) {
        std::cerr << "AdjugateDeterminant: non-exact back-substitution\n";
        throw TerminalException{1};
      }
#endif
      xhat[i] = quot;
    }
    for (int i = 0; i < n; i++)
      Output(i, c) = xhat[i];
  }
  return {std::move(Output), std::move(det)};
}

// Matrix inverse through the fraction-free LU factorization:
// A^{-1} = adj(A) / det(A). Forward-then-back does fewer ring operations
// than the full Gauss-Jordan, so it is the cheaper fraction-free inverse.
template <typename T> MyMatrix<T> InverseFractionFreeLU(MyMatrix<T> const &Input) {
  std::pair<MyMatrix<T>, T> pair = AdjugateDeterminant(Input);
  MyMatrix<T> const &adj = pair.first;
  T const &det = pair.second;
  int n = Input.rows();
  MyMatrix<T> Output(n, n);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      T quot = adj(i, j) / det;
#ifdef SANITY_CHECK_MATRIX_INVERSE
      if (quot * det != adj(i, j)) {
        std::cerr << "InverseFractionFreeLU: A^{-1} is not representable over T "
                     "(the determinant does not divide the adjugate)\n";
        throw TerminalException{1};
      }
#endif
      Output(i, j) = quot;
    }
  return Output;
}

// Types opting into the fraction-free LU inverse (see use_fraction_free_lu): the
// integer rings and the exact number fields. Forward-elimination-then-back-
// substitution does fewer ring operations than the classical Gauss-Jordan.
template <typename T>
requires (use_fraction_free_lu<T>::value)
inline MyMatrix<T> Inverse(MyMatrix<T> const &Input) {
  return InverseFractionFreeLU(Input);
}

// Ordinary field not opting into fraction-free LU (mpq_class and the rational
// fields, floating point, Fp, ...): classical Gauss-Jordan, whose
// SelectBestPivot pivoting also preserves numerical stability for floats.
template <typename T>
requires (!use_fraction_free_lu<T>::value && is_ring_field<T>::value)
inline MyMatrix<T> Inverse(MyMatrix<T> const &Input) {
  return InverseKernel(Input);
}

// Any other non-field ring not opting into fraction-free LU: map to the
// overlying field and invert there.
template <typename T>
requires (!use_fraction_free_lu<T>::value && !is_ring_field<T>::value)
inline MyMatrix<T> Inverse(MyMatrix<T> const &Input) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> InputF = UniversalMatrixConversion<Tfield, T>(Input);
  MyMatrix<Tfield> OutputF = InverseKernel(InputF);
  return UniversalMatrixConversion<T, Tfield>(OutputF);
}

template <typename T> MyMatrix<T> CongrMap(MyMatrix<T> const &eMat) {
  MyMatrix<T> TheInv = Inverse(eMat);
  return TransposedMat(TheInv);
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXINVERSE_H_
// clang-format on
