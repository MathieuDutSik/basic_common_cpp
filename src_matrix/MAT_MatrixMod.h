#ifndef INCLUDE_MAT_MATRIX_MOD_H
#define INCLUDE_MAT_MATRIX_MOD_H

#include "MAT_MatrixInt.h"


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
MyMatrix<T> ComputeBasisInvariantSpace(std::vector<MyMatrix<T>> const& ListMat, MyMatrix<T> const& TheSpace, T const& TheMod)
{
  int n = TheSpace.rows();
  if (ListMat.size() == 0)
    return TheSpace;
  int n_mat = ListMat.size();
  MyMatrix<T> Equa(2 * n, n_mat * n);
  for (int i_mat=0; i_mat<n_mat; i_mat++) {
    MyMatrix<T> const& eMat = ListMat[i_mat];
    MyMatrix<T> eProd = TheSpace * eMat - TheSpace;
    for (int i=0; i<n; i++) {
      for (int j=0; j<n; j++) {
        Equa(i,j + i_mat * n) = eProd(i,j);
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
  MyMatrix<T> NSP_red(n,n);
  for (int i=0; i<n_row; i++)
    for (int j=0; j<n; j++)
      NSP_red(i,j) = NSP(i,j);
  return NSP_red * TheSpace;
}

/*
Equation to solve is x1 M1 = x2 M2 + u MOD
*/
template <typename T>
MyMatrix<T> IntersectionLatticeMod(MyMatrix<T> const& M1, MyMatrix<T> const& M2, T const& TheMod)
{
  int n_row1 = M1.rows();
  int n_row2 = M2.rows();
  int n = M1.cols();
  MyMatrix<T> Equa(n_row1 + n_row2 + n, n);
  for (int i1=0; i1<n_row1; i1++)
    for (int j=0; j<n; j++)
      Equa(i1, j) = M1(i1, j);
  for (int i2=0; i2<n_row2; i2++)
    for (int j=0; j<n; j++)
      Equa(n_row1 + i2, j) = M1(i2, j);
  for (int i=0; i<n; i++)
    for (int j=0; j<n; j++) {
      if (i == j) {
        Equa(n_row1 + n_row2 + i, j) = TheMod;
      } else {
        Equa(n_row1 + n_row2 + i, j) = 0;
      }
    }
  MyMatrix<T> NSP = NullspaceIntMat(Equa);
  MyMatrix<T> NSP_red(NSP.rows(), n_row1);
  for (int i=0; i<NSP.rows(); i++)
    for (int j=0; j<n_row1; j++)
      NSP_red(i,j) = NSP(i,j);
  return NSP_red * M1;
}




#endif
