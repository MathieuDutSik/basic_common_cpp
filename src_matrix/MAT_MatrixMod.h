// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXMOD_H_
#define SRC_MATRIX_MAT_MATRIXMOD_H_

// clang-format off
#include "MAT_MatrixInt.h"
#include "quadratic_residue.h"
#include "factorizations.h"
#include "NumberTheoryGeneric.h"
#include <vector>
#include <cmath>
// clang-format on

template <typename T, typename Tmod>
MyMatrix<Tmod> ModuloReductionMatrix(MyMatrix<T> const &M, T const &TheMod) {
  int n_row = M.rows();
  int n_col = M.cols();
  MyMatrix<Tmod> RetMat(n_row, n_col);
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++) {
      T val = ResInt(M(i, j), TheMod);
      RetMat(i, j) = UniversalScalarConversion<Tmod, T>(val);
    }
  }
  return RetMat;
}

template <typename T, typename Tmod>
std::vector<MyMatrix<Tmod>>
ModuloReductionStdVectorMatrix(std::vector<MyMatrix<T>> const &ListM,
                               T const &TheMod) {
  std::vector<MyMatrix<Tmod>> ListRetMat;
  for (auto &M : ListM)
    ListRetMat.push_back(ModuloReductionMatrix<T, Tmod>(M, TheMod));
  return ListRetMat;
}

template <typename T, typename Tmod>
MyVector<Tmod> ModuloReductionVector(MyVector<T> const &V, T const &TheMod) {
  int siz = V.size();
  MyVector<Tmod> retV(siz);
  for (int i = 0; i < siz; i++) {
    T val = ResInt(V(i), TheMod);
    retV(i) = UniversalScalarConversion<Tmod, T>(val);
  }
  return retV;
}

template <typename T>
MyVector<T> VectorMod(MyVector<T> const &V, T const &TheMod) {
  int n = V.size();
  MyVector<T> Vret(n);
  for (int i = 0; i < n; i++)
    Vret(i) = ResInt(V(i), TheMod);
  return Vret;
}

/*
We want to find the vectors x in Z^n such that
x TheSpace P = x TheSpace + u MOD
 */
template <typename T>
MyMatrix<T> ComputeBasisInvariantSpace(std::vector<MyMatrix<T>> const &ListMat,
                                       MyMatrix<T> const &TheSpace,
                                       T const &TheMod) {
  int n = TheSpace.rows();
  if (ListMat.size() == 0)
    return TheSpace;
  int n_mat = ListMat.size();
  MyMatrix<T> Equa(2 * n, n_mat * n);
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    MyMatrix<T> const &eMat = ListMat[i_mat];
    MyMatrix<T> eProd = TheSpace * eMat - TheSpace;
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        Equa(i, j + i_mat * n) = eProd(i, j);
        if (i == j) {
          Equa(i + n, j + i_mat * n) = TheMod;
        } else {
          Equa(i + n, j + i_mat * n) = 0;
        }
      }
    }
  }
  MyMatrix<T> NSP = NullspaceIntMat(Equa);
  int n_row = NSP.rows();
  MyMatrix<T> NSP_red(n, n);
  for (int i = 0; i < n_row; i++)
    for (int j = 0; j < n; j++)
      NSP_red(i, j) = NSP(i, j);
  return NSP_red * TheSpace;
}

/*
Equation to solve is x1 M1 = x2 M2 + u MOD
*/
template <typename T>
MyMatrix<T> IntersectionLatticeMod(MyMatrix<T> const &M1, MyMatrix<T> const &M2,
                                   T const &TheMod) {
  int n_row1 = M1.rows();
  int n_row2 = M2.rows();
  int n = M1.cols();
  MyMatrix<T> Equa(n_row1 + n_row2 + n, n);
  for (int i1 = 0; i1 < n_row1; i1++)
    for (int j = 0; j < n; j++)
      Equa(i1, j) = M1(i1, j);
  for (int i2 = 0; i2 < n_row2; i2++)
    for (int j = 0; j < n; j++)
      Equa(n_row1 + i2, j) = M1(i2, j);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++) {
      if (i == j) {
        Equa(n_row1 + n_row2 + i, j) = TheMod;
      } else {
        Equa(n_row1 + n_row2 + i, j) = 0;
      }
    }
  MyMatrix<T> NSP = NullspaceIntMat(Equa);
  MyMatrix<T> NSP_red(NSP.rows(), n_row1);
  for (int i = 0; i < NSP.rows(); i++)
    for (int j = 0; j < n_row1; j++)
      NSP_red(i, j) = NSP(i, j);
  return NSP_red * M1;
}

template <typename T>
bool Kernel_FindIsotropicVectorModQuadResidue(MyMatrix<T> const &M,
                                              MyVector<T> &V, T const &TheMod) {
  int n = M.rows();
  T cst(0);
  for (int i = 0; i < n - 1; i++) {
    for (int j = 0; j < n - 1; j++) {
      cst += M(i, j) * V(i) * V(j);
    }
  }
  T lin(0);
  for (int i = 0; i < n - 1; i++) {
    lin += M(i, n - 1) * V(i);
  }
  T C = M(n - 1, n - 1);
  // The equation to be solved becomes
  // 0 = cst + 2 lin x_n + C x_n^2
  // If C = 0 then solution is obvious.
  // Otherwise, divide by C and get
  // 0 = cst + 2 lin x_n + x_n^2
  // 0 = cst + (x_n + lin)^2 - lin^2
  // lin^2 - cst = (x_n + lin)^2
  T Cinv = mod_inv(C, TheMod);
  lin *= Cinv;
  cst *= Cinv;
  T a = lin * lin - cst;
  std::optional<T> opt = find_quadratic_residue(a, TheMod);
  if (opt) {
    T const &y = *opt;
    T xn = y - lin;
    V(n - 1) = xn;
#ifdef DEBUG_MOD_OPERATION
    T sum = 0;
    for (int i = 0; i < n; i++)
      for (int j = 0; j < n; j++)
        sum += M(i, j) * V(i) * V(j);
    T res = ResInt(sum, TheMod);
    if (res != 0) {
      std::cerr << "We failed to find an isotropic vector\n";
      throw TerminalException{1};
    }
#endif
    return true;
  }
  return false;
}

/*
  Find isotropic vector
  We follow here the paper "Quadratic equations in dimensions 4, 5 and more"
  and in particular Lemma 1.
*/
template <typename T>
MyVector<T> FindIsotropicVectorModRandom(MyMatrix<T> const &M,
                                         T const &TheMod) {
  int n = M.rows();
  MyVector<T> V = ZeroVector<T>(n);
  for (int i = 0; i < n; i++) {
    if (M(i, i) == 0) {
      V(i) = 1;
      return V;
    }
  }
  while (true) {
    // We set up the first n-1 coordinates at random
    // Then we solve the equation for finding the last one.
    for (int i = 0; i < n - 1; i++) {
      int val = rand();
      T val_T(val);
      V(i) = val_T;
    }
    bool test = Kernel_FindIsotropicVectorModQuadResidue(M, V, TheMod);
    if (test) {
      return V;
    }
  }
}

/*
  About the solution for two-dimensional.
 */
template <typename T>
std::optional<MyVector<T>> FindIsotropicVectorModTwoDim(MyMatrix<T> const &M,
                                                        T const &TheMod) {
  MyVector<T> V(2);
  if (M(1, 1) == 0) {
    V(0) = 0;
    V(1) = 1;
    return V;
  }
  V(0) = 1;
  bool test = Kernel_FindIsotropicVectorModQuadResidue(M, V, TheMod);
  if (test) {
    return V;
  }
  return {};
}

template <typename T>
std::optional<MyVector<T>> FindIsotropicVectorMod_Z(MyMatrix<T> const &M,
                                                    T const &TheMod) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  int n = M.rows();
  MyMatrix<T> Mred(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Mred(i, j) = ResInt(M(i, j), TheMod);
    }
  }
  if (n == 1) {
    T val = Mred(0, 0);
    if (val == 0) {
      MyVector<T> V(1);
      V(0) = 1;
      return V;
    }
    return {};
  }
  if (n == 2) {
    return FindIsotropicVectorModTwoDim(Mred, TheMod);
  }
  return FindIsotropicVectorModRandom(Mred, TheMod);
}

template <typename T>
inline typename std::enable_if<is_implementation_of_Q<T>::value,
                               std::optional<MyVector<T>>>::type
FindIsotropicVectorMod(MyMatrix<T> const &M, T const &TheMod) {
  using Tring = typename underlying_ring<T>::ring_type;
  MyMatrix<T> M1 = RemoveFractionMatrix(M);
  MyMatrix<Tring> M2 = UniversalMatrixConversion<Tring, T>(M1);
  Tring TheMod_ring = UniversalScalarConversion<Tring, T>(TheMod);
  std::optional<MyVector<Tring>> opt =
      FindIsotropicVectorMod_Z(M2, TheMod_ring);
  if (opt) {
    MyVector<Tring> const &eV = *opt;
    MyVector<T> eV_T = UniversalVectorConversion<T, Tring>(eV);
    return eV_T;
  }
  return {};
}

template <typename T>
inline typename std::enable_if<!is_implementation_of_Q<T>::value,
                               std::optional<MyVector<T>>>::type
FindIsotropicVectorMod(MyMatrix<T> const &M, T const &TheMod) {
  return FindIsotropicVectorMod_Z(M, TheMod);
}

template <typename T>
bool IsMatrixZeroMod(MyMatrix<T> const &M, T const &TheMod) {
  int n_row = M.rows();
  int n_col = M.cols();
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++) {
      T res = ResInt(M(i, j), TheMod);
      if (res != 0) {
        return false;
      }
    }
  }
  return true;
}

/*
  We want to find the solution of M x = 0
  with the operation being modulo p.

  We tried before a number of tricks by extending the matrix
  but that seems not to work.
 */
template <typename T>
SelectionRowCol<T> SelectRowColMatMod(MyMatrix<T> const &M, T const &TheMod) {
  int n_row = M.rows();
  int n_col = M.cols();
  int maxRank = n_row;
  if (n_col < maxRank)
    maxRank = n_col;
  size_t sizMat = maxRank + 1;
  MyMatrix<T> Mwork(sizMat, n_col);
  int miss_val = -1;
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++) {
      Mwork(i, j) = ResInt(M(i, j), TheMod);
    }
  }
  std::vector<int> ListColSelect;
  std::vector<uint8_t> ListColSelect01(n_col, 0);
  std::vector<int> ListRowSelect;
  int eRank = 0;
  for (int iRow = 0; iRow < n_row; iRow++) {
#ifdef DEBUG_MATRIX_MOD
    std::cerr << "MATMOD: NullspaceTrMatMod iRow=" << iRow << "\n";
#endif
    for (int iCol = 0; iCol < n_col; iCol++) {
      T res = ResInt(M(iRow, iCol), TheMod);
      Mwork(eRank, iCol) = res;
    }
    for (int iRank = 0; iRank < eRank; iRank++) {
      int eCol = ListColSelect[iRank];
      T eVal1 = Mwork(eRank, eCol);
      if (eVal1 != 0) {
        for (int iCol = 0; iCol < n_col; iCol++) {
          T val = Mwork(eRank, iCol) - eVal1 * Mwork(iRank, iCol);
          Mwork(eRank, iCol) = ResInt(val, TheMod);
        }
      }
    }
    auto get_firstnonzero_iife = [&]() -> int {
      for (int iCol = 0; iCol < n_col; iCol++) {
        if (Mwork(eRank, iCol) != 0) {
          return iCol;
        }
      }
      return miss_val;
    };
    int FirstNZ = get_firstnonzero_iife();
    if (FirstNZ != miss_val) {
      ListColSelect.push_back(FirstNZ);
      ListColSelect01[FirstNZ] = 1;
      ListRowSelect.push_back(iRow);
      T eVal1 = Mwork(eRank, FirstNZ);
      T eVal2 = mod_inv(eVal1, TheMod);
#ifdef DEBUG_MATRIX_MOD
      std::cerr << "MATMOD: eVal1=" << eVal1 << " eVal2=" << eVal2 << "\n";
#endif
      for (int iCol = 0; iCol < n_col; iCol++) {
        T val = Mwork(eRank, iCol) * eVal2;
        Mwork(eRank, iCol) = ResInt(val, TheMod);
      }
      for (int iRank = 0; iRank < eRank; iRank++) {
        T eVal1 = Mwork(iRank, FirstNZ);
        if (eVal1 != 0) {
          int StartCol = ListColSelect[iRank];
          for (int iCol = StartCol; iCol < n_col; iCol++) {
            T val = Mwork(iRank, iCol) - eVal1 * Mwork(eRank, iCol);
            Mwork(iRank, iCol) = ResInt(val, TheMod);
          }
        }
      }
      eRank++;
    }
  }
#ifdef DEBUG_MATRIX_MOD
  std::cerr << "MATMOD: After the main loop\n";
#endif
  int nbVectZero = n_col - eRank;
#ifdef DEBUG_MATRIX_MOD
  std::cerr << "MATMOD: n_col=" << n_col << " eRank=" << eRank
            << " nbVectZero=" << nbVectZero << "\n";
  std::cerr << "MATMOD: ListColSelect01=";
  for (auto &k : ListColSelect01) {
    std::cerr << static_cast<int>(k);
  }
  std::cerr << "\n";
#endif
  MyMatrix<T> NSP = ZeroMatrix<T>(nbVectZero, n_col);
  int nbVect = 0;
  for (int iCol = 0; iCol < n_col; iCol++) {
    if (ListColSelect01[iCol] == 0) {
      NSP(nbVect, iCol) = TheMod - 1;
      for (int iRank = 0; iRank < eRank; iRank++) {
        int eCol = ListColSelect[iRank];
        NSP(nbVect, eCol) = Mwork(iRank, iCol);
      }
      nbVect += 1;
    }
  }
#ifdef DEBUG_MATRIX_MOD
  std::cerr << "MATMOD: NSP built\n";
  MyMatrix<T> prod = M * NSP.transpose();
  if (!IsMatrixZeroMod(prod, TheMod)) {
    std::cerr << "SelectRowColMatMod: The matrix prod is not zero\n";
    throw TerminalException{1};
  }
#endif
  return {static_cast<size_t>(eRank), NSP, ListColSelect, ListRowSelect};
}

template <typename T>
MyMatrix<T> NullspaceTrMatMod(MyMatrix<T> const &M, T const &TheMod) {
  SelectionRowCol<T> result = SelectRowColMatMod(M, TheMod);
  return result.NSP;
}

template <typename T>
MyMatrix<T> NullspaceMatMod(MyMatrix<T> const &M, T const &TheMod) {
  MyMatrix<T> Mtr = M.transpose();
  MyMatrix<T> NSP = NullspaceTrMatMod(Mtr, TheMod);
#ifdef DEBUG_MATRIX_MOD
  MyMatrix<T> prod = NSP * M;
  if (!IsMatrixZeroMod(prod, TheMod)) {
    std::cerr << "NullspaceMatProd: The matrix prod is not zero\n";
    throw TerminalException{1};
  }
#endif
  return NSP;
}

template <typename T>
T DeterminantMatMod(MyMatrix<T> const &TheMat, T const &TheMod) {
  int n = TheMat.rows();
  if (n != TheMat.cols()) {
    std::cerr << "DeterminantMatMod: Matrix must be square\n";
    throw TerminalException{1};
  }

  if (n == 0) {
    return T(1);
  }

  if (n == 1) {
    return ResInt(TheMat(0, 0), TheMod);
  }

  MyMatrix<T> M(n, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      M(i, j) = ResInt(TheMat(i, j), TheMod);
    }
  }

  T det = T(1);
  for (int i = 0; i < n; i++) {
    // Find pivot
    int pivot_row = -1;
    for (int k = i; k < n; k++) {
      if (M(k, i) != 0) {
        pivot_row = k;
        break;
      }
    }

    if (pivot_row == -1) {
      // No pivot found, determinant is zero
      return T(0);
    }

    // Swap rows if needed
    if (pivot_row != i) {
      for (int j = 0; j < n; j++) {
        T temp = M(i, j);
        M(i, j) = M(pivot_row, j);
        M(pivot_row, j) = temp;
      }
      det = TheMod - det;
    }

    T pivot = M(i, i);
    det *= pivot;
    det = ResInt(det, TheMod);

    // Get multiplicative inverse of pivot
    T pivot_inv = mod_inv(pivot, TheMod);

    // Eliminate column
    for (int k = i + 1; k < n; k++) {
      if (M(k, i) != 0) {
        T val = M(k, i) * pivot_inv;
        T factor = ResInt(val, TheMod);
        for (int j = i; j < n; j++) {
          T val1 = factor * M(i, j);
          T val = M(k, j) - ResInt(val1, TheMod);
          M(k, j) = ResInt(val, TheMod);
        }
      }
    }
  }

  return det;
}

template <typename T>
MyMatrix<T> SmithNormalFormIntegerMat(MyMatrix<T> const &TheMat) {
  static_assert(is_ring_field<T>::value || is_euclidean_domain<T>::value,
                "Requires T to be a ring or euclidean domain");

  int n = TheMat.rows();
  int m = TheMat.cols();

  if (n == 0 || m == 0) {
    return TheMat;
  }

  MyMatrix<T> A = TheMat;
  int rank = 0;
  int minDim = std::min(n, m);

  for (int k = 0; k < minDim; k++) {
    // Find pivot with smallest non-zero norm
    int pivot_row = -1;
    int pivot_col = -1;
    T min_norm = T(0);
    bool found_pivot = false;

    for (int i = k; i < n; i++) {
      for (int j = k; j < m; j++) {
        if (A(i, j) != 0) {
          T norm = T_abs(A(i, j));
          if (!found_pivot || norm < min_norm) {
            min_norm = norm;
            pivot_row = i;
            pivot_col = j;
            found_pivot = true;
          }
        }
      }
    }

    if (!found_pivot) {
      break; // No more pivots
    }

    // Move pivot to position (k, k)
    if (pivot_row != k) {
      for (int j = 0; j < m; j++) {
        T temp = A(k, j);
        A(k, j) = A(pivot_row, j);
        A(pivot_row, j) = temp;
      }
    }
    if (pivot_col != k) {
      for (int i = 0; i < n; i++) {
        T temp = A(i, k);
        A(i, k) = A(i, pivot_col);
        A(i, pivot_col) = temp;
      }
    }

    // Make pivot positive
    if (A(k, k) < 0) {
      for (int j = k; j < m; j++) {
        A(k, j) = -A(k, j);
      }
    }

    // Clear column k below diagonal
    for (int i = k + 1; i < n; i++) {
      if (A(i, k) != 0) {
        T g = GcdPair(A(k, k), A(i, k));
        T u = A(k, k) / g;
        T v = A(i, k) / g;

        for (int j = k; j < m; j++) {
          T temp = u * A(i, j) - v * A(k, j);
          A(k, j) = A(k, j) + A(i, j);
          A(i, j) = temp;
        }
      }
    }

    // Clear row k to the right of diagonal
    for (int j = k + 1; j < m; j++) {
      if (A(k, j) != 0) {
        T g = GcdPair(A(k, k), A(k, j));
        T u = A(k, k) / g;
        T v = A(k, j) / g;

        for (int i = k; i < n; i++) {
          T temp = u * A(i, j) - v * A(i, k);
          A(i, k) = A(i, k) + A(i, j);
          A(i, j) = temp;
        }
      }
    }

    rank++;

    // Check if we need to continue (if off-diagonal elements were created)
    bool need_continue = false;
    for (int i = k + 1; i < n && !need_continue; i++) {
      if (A(i, k) != 0) need_continue = true;
    }
    for (int j = k + 1; j < m && !need_continue; j++) {
      if (A(k, j) != 0) need_continue = true;
    }

    if (need_continue) {
      k--; // Redo this step
    }
  }

  // Ensure divisibility condition: d[i] divides d[i+1]
  for (int i = 0; i < rank - 1; i++) {
    if (A(i, i) != 0 && A(i + 1, i + 1) != 0) {
      T g = GcdPair(A(i, i), A(i + 1, i + 1));
      if (g != A(i, i)) {
        // Need to adjust to maintain divisibility
        A(i, i) = g;
        A(i + 1, i + 1) = (A(i, i) * A(i + 1, i + 1)) / g;
      }
    }
  }

  return A;
}

template <typename T>
double HadamardUpperBound(MyMatrix<T> const &TheMat) {
  int n = TheMat.rows();
  if (n != TheMat.cols()) {
    std::cerr << "HadamardUpperBound: Matrix must be square\n";
    throw TerminalException{1};
  }
  
  if (n == 0) {
    return 1.0;
  }
  
  if (n == 1) {
    return UniversalScalarConversion<double, T>(T_abs(TheMat(0, 0)));
  }
  
  double bound = 1.0;
  
  for (int i = 0; i < n; i++) {
    double row_norm_squared = 0.0;
    
    // Compute squared Euclidean norm of row i
    for (int j = 0; j < n; j++) {
      double val = UniversalScalarConversion<double, T>(TheMat(i, j));
      row_norm_squared += val * val;
    }
    
    // Take square root to get the actual norm
    double row_norm = std::sqrt(row_norm_squared);
    
    // Multiply to the bound
    bound *= row_norm;
  }
  
  return bound;
}

template <typename T>
T DeterminantMatHadamard(MyMatrix<T> const &TheMat) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  
  int n = TheMat.rows();
  if (n != TheMat.cols()) {
    std::cerr << "DeterminantMatHadamard: Matrix must be square\n";
    throw TerminalException{1};
  }
  
  if (n == 0) {
    return T(1);
  }
  
  if (n == 1) {
    return TheMat(0, 0);
  }
  
  // Step 1: Compute Hadamard upper bound
  double bound = HadamardUpperBound(TheMat);
  T target_product = T(3) * T(static_cast<long long>(std::ceil(bound)));
  
#ifdef DEBUG_MATRIX_MOD
  std::cerr << "DETHADAMARD: Hadamard bound = " << bound << "\n";
  std::cerr << "DETHADAMARD: Target product = " << target_product << "\n";
#endif
  
  // Step 2: Generate primes and compute determinants modulo primes
  PrimeGenerator<T> prime_gen;
  std::vector<T> primes;
  std::vector<T> det_mods;
  T product(1);
  
  while (product < target_product) {
    T prime = prime_gen.get_prime();
    primes.push_back(prime);
    
    T det_mod = DeterminantMatMod(TheMat, prime);
    det_mods.push_back(det_mod);
    
    product *= prime;
    
#ifdef DEBUG_MATRIX_MOD
    std::cerr << "DETHADAMARD: Prime = " << prime << ", det mod = " << det_mod << ", product = " << product << "\n";
#endif
  }
  
#ifdef DEBUG_MATRIX_MOD
  std::cerr << "DETHADAMARD: Used " << primes.size() << " primes\n";
  std::cerr << "DETHADAMARD: Final product = " << product << "\n";
#endif
  
  // Step 3: Use Chinese Remainder Theorem to find determinant mod product
  T det_crt = chinese_remainder_theorem(det_mods, primes);
  
  // Step 4: Find the value nearest to 0
  T half_product = product / T(2);
  
  // If det_crt > product/2, then det_crt - product is closer to 0
  if (det_crt > half_product) {
    det_crt = det_crt - product;
  }
  
#ifdef DEBUG_MATRIX_MOD
  std::cerr << "DETHADAMARD: Final determinant = " << det_crt << "\n";
#endif
  
  return det_crt;
}

template <typename T>
int DimensionKernelMod(MyMatrix<T> const &M, T const &TheMod) {
  int n_row = M.rows();
  int n_col = M.cols();
  MyMatrix<T> Mbig(n_row + n_col, n_col);
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++) {
      Mbig(i, j) = M(i, j);
    }
  }
  for (int i = 0; i < n_col; i++) {
    for (int j = 0; j < n_col; j++) {
      T val(0);
      if (i == j)
        val = TheMod;
      Mbig(i + n_row, j) = val;
    }
  }
  MyMatrix<T> NSP = NullspaceIntMat(Mbig);
  return NSP.rows();
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXMOD_H_
// clang-format on
