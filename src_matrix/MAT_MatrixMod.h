// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXMOD_H_
#define SRC_MATRIX_MAT_MATRIXMOD_H_

// clang-format off
#include "MAT_MatrixInt.h"
#include "quadratic_residue.h"
#include <vector>
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

template <typename T> T mod_inv(T const &a, T const &P) {
  T t = 0;
  T newt = 1;
  T r = P;
  T newr = m;
  T q, tmp;
  while (newr != 0) {
    q = r / newr;
    tmp = t;
    t = newt;
    newt = tmp - q * newt;
    tmp = r;
    r = newr;
    newr = tmp - q * newr;
  }
  if (r > 1)
    return 0;
  if (t < 0)
    t = t + P;
  return t;
}


template<typename T>
bool Kernel_FindIsotropicVectorQuadResidue(MyMatrix<T> const &M, MyVector<T> & V, T const &TheMod) {
  int n = M.rows();
  T cst = 0;
  for (int i = 0; i < n - 1; i++) {
    for (int j = 0; j < n - 1; j++) {
      cst += M(i, j) * V(i) * V(j);
    }
  }
  T lin = 0;
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
      throw TerminalExcpetion{1};
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
template<typename T>
MyVector<T> FindIsotropicVectorRandom(MyMatrix<T> const &M, T const &TheMod) {
  int n = M.rows();
  MyVector<T> V(n);
  if (M(n-1,n-1) == 0) {
    for (int i=0; i<n-1; i++)
      V(i) = 0;
    V(n-1) = 1;
    return true;
  }
  while (true) {
    // We set up the first n-1 coordinates at random
    // Then we solve the equation for finding the last one.
    for (int i = 0; i < n - 1; i++) {
      int val = rand();
      T val_T(val);
      V(i) = val_T;
    }
    bool test = Kernel_FindIsotropicVectorQuadResidue(M, V, TheMod);
    if (test) {
      return V;
    }
  }
}

/*
  About the solution for two-dimensional.
 */
template<typename T>
std::optional<MyVector<T>> FindIsotropicVectorTwoDim(MyMatrix<T> const &M, T const &TheMod) {
  MyVector<T> V(2);
  if (M(1,1) == 0) {
    V(0) = 0;
    V(1) = 1;
    return V;
  }
  V(0) = 1;
  bool test = Kernel_FindIsotropicVectorQuadResidue(M, V, TheMod);
  if (test) {
    return V;
  }
  return {};
}

template<typename T>
std::optional<MyVector<T>> FindIsotropicVector(MyMatrix<T> const &M, T const &TheMod) {
  int n = M.rows();
  if (n == 1) {
    T val = M(0,0);
    if (val == 0) {
      MyVector<T> V(1);
      V(0) = 1;
      return V;
    }
    return {};
  }
  if (n == 2) {
    return FindIsotropicVectorTwoDim(M, TheMod);
  }
  return FindIsotropicVectorRandom(M, TheMod);
}





template<typename T>
struct ResultNullspaceMod {
  int dimNSP;
  MyMatrix<T> BasisTot;
};


/*
  We work with the matrix M with m rows and n cols.
  We want to solve the equation x M = 0 modulo p. So, x has m cols.
  We can rewrite the equation over Z as:
  xM + y p I_n = 0 with y a vector with n columns.
  This gets us
  (x, y) / M     \
         \ p I_n /
  And so we can find a basis of the kernel.

  What we want to have is a Z-basis (x^1, ...., x^d, .... x^m)
  with (x^1, ...., x^d) image in (F_p)^m is spanning the
  kernel over Fp of M.


 */
template<typename T>
ResultNullspaceMod<T> NullspaceMatMod(MyMatrix<T> const& M, T const& TheMod) {
  int n_row = M.rows();
  int n_col = M.cols();
  MyMatrix<T> Mbig(n_row + n_col, n_col);
  for (int i=0; i<n_row; i++) {
    for (int j=0; j<n_col; j++) {
      Mbig(i, j) = M(i, j);
    }
  }
  for (int i=0; i<n_col; i++) {
    for (int j=0; j<n_col; j++) {
      T val(0);
      if (i == j)
        val = TheMod;
      Mbig(i + n_row, j) = val;
    }
  }
  // 1) compute the full null-space
  MyMatrix<T> NSP = NullspaceIntMat(Mbig);
  int dimNSP = NSP.rows();
  if (dimNSP == n_row) {
    return IdentityMat<T>(n_row);
  }
  if (dimNSP == 0) {
    MyMatrix<T> Mret(0, n_row);
    return Mret;
  }
  // 2) Compute the reduced null-space (but it could
  MyMatrix<T> NSPred(dimNSP, n_row);
  for (int i=0; i<dimNSP; i++) {
    for (int j=0; j<n_row; j++) {
      NSPred(i, j) = NSP(i, j);
    }
  }
#ifdef DEBUG_MATRIX_MOD
  dim rnkA = RankMat(NSPred);
  if (rnkA != dimNSP) {
    std::cerr << "We have rnkA=" << rnkA << " dimNSP=" << dimNSP << "\n";
    std::cerr << "This is not what we expected\n";
    throw TerminalException{1};
  }
  MyMatrix<T> prod = NSPred * M;
  for (int i=0; i<prod.rows(); i++) {
    for (int j=0; j<prod.cols(); j++) {
      T res = ResInt(prod(i,j), TheMod);
      if (res != 0) {
        std::cerr << "The residue is not 0\n";
        throw TerminalException{1};
      }
    }
  }
#endif
  MyMatrix<T> BasisComp = SubspaceCompletionInt(NSPred, n_row);
  MyMatrix<T> FullMat = Concatenate(NSPred, BasisComp);
  return {dimNSP, std::move(FullMat)};
}

template<typename T>
int DimensionKernelMod(MyMatrix<T> const& M, T const& TheMod) {
  int n_row = M.rows();
  int n_col = M.cols();
  MyMatrix<T> Mbig(n_row + n_col, n_col);
  for (int i=0; i<n_row; i++) {
    for (int j=0; j<n_col; j++) {
      Mbig(i, j) = M(i, j);
    }
  }
  for (int i=0; i<n_col; i++) {
    for (int j=0; j<n_col; j++) {
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
