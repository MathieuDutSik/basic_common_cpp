// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIXSOLUTIONMAT_H_
#define SRC_MATRIX_MAT_MATRIXSOLUTIONMAT_H_

// Solving the linear system Y = X A and the use cases built on it.
//
//  --- SolutionMatKernel / SolutionMat: one solution X of Y = X A if it exists.
//  --- SolutionMatRepetitive: shares the preprocessing across many right-hand
//      sides.
//  --- TestEqualitySpannedSpaces, ListSolutionMat,
//      ExpressVectorsInIndependentFamilt: utilities expressed via SolutionMat.
//  --- SolutionMat_LeastSquare: floating-point least-squares solve (via Eigen).

// clang-format off
#include "MAT_MatrixInverse.h"
#include "MAT_MatrixRankmat.h"
#include <optional>
#include <vector>
// clang-format on

// Given the equation Y = XA, we find one solution X if it exists.
//
template <typename T>
std::optional<MyVector<T>> SolutionMatKernel(MyMatrix<T> const &eMat,
                                             MyVector<T> const &eVect) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in SolutionMat");
  if (eMat.rows() == 0) {
    if (!IsZeroVector(eVect))
      return {};
    MyVector<T> eSol(0);
    return eSol;
  }
  int nbRow = eMat.rows();
  int nbCol = eMat.cols();
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(eMat);
  int eRank = eSelect.TheRank;
  std::vector<int> ListRowSelect = eSelect.ListRowSelect;
  std::vector<int> ListColSelect = eSelect.ListColSelect;
  MyMatrix<T> eMat2 = SelectRow(eMat, ListRowSelect);
  MyMatrix<T> eMat3 = SelectColumn(eMat2, ListColSelect);
  MyVector<T> eVectB = SelectColumnVector(eVect, ListColSelect);
  MyMatrix<T> eMatInv = Inverse(eMat3);
  MyVector<T> eSol = ProductVectorMatrix(eVectB, eMatInv);
  MyVector<T> eProd = ProductVectorMatrix(eSol, eMat2);
  for (int iCol = 0; iCol < nbCol; iCol++)
    if (eProd(iCol) != eVect(iCol))
      return {};
  MyVector<T> eRetSol = ZeroVector<T>(nbRow);
  for (int iRank = 0; iRank < eRank; iRank++) {
    int iRow = ListRowSelect[iRank];
    eRetSol(iRow) = eSol(iRank);
  }
  return eRetSol;
}

template <typename T>
requires is_ring_field<T>::value
inline std::optional<MyVector<T>> SolutionMat(MyMatrix<T> const &eMat,
                                              MyVector<T> const &eVect) {
  return SolutionMatKernel(eMat, eVect);
}

// Over a ring the solution is fractional in general, so it is returned
// over the overlying field; the return is empty exactly when the system
// has no solution at all. For a solution over the ring itself, use
// SolutionIntMat.
template <typename T>
requires (!is_ring_field<T>::value)
inline std::optional<MyVector<typename overlying_field<T>::field_type>>
SolutionMat(MyMatrix<T> const &eMat, MyVector<T> const &eVect) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> eMatF = UniversalMatrixConversion<Tfield, T>(eMat);
  MyVector<Tfield> eVectF = UniversalVectorConversion<Tfield, T>(eVect);
  return SolutionMatKernel(eMatF, eVectF);
}

template <typename T> struct SolutionMatRepetitive {
private:
  MyMatrix<T> TheBasis;
  std::vector<int> ListColSelect;
  MyMatrix<T> InvMat;

public:
  SolutionMatRepetitive(MyMatrix<T> const &_TheBasis) : TheBasis(_TheBasis) {
#ifdef SANITY_CHECK_MAT_MATRIX
    if (RankMat(TheBasis) != TheBasis.rows()) {
      std::cerr << "RankMat(TheBasis)=" << RankMat(TheBasis) << "\n";
      std::cerr << "  TheBasis.rows()=" << TheBasis.rows() << "\n";
      std::cerr << "Error in SolutionMatRepetitive\n";
      throw TerminalException{1};
    }
#endif
    SelectionRowCol<T> src = TMat_SelectRowCol(TheBasis);
    ListColSelect = src.ListColSelect;
    MyMatrix<T> SelMat = SelectColumn(TheBasis, ListColSelect);
    InvMat = Inverse(SelMat);
  }
  std::optional<MyVector<T>> GetSolution(MyVector<T> const &eVect) {
    int siz = ListColSelect.size();
    MyVector<T> V(siz);
    for (int u = 0; u < siz; u++) {
      V(u) = eVect(ListColSelect[u]);
    }
    MyVector<T> MySol = InvMat.transpose() * V;
    if (TheBasis.transpose() * MySol != eVect) {
      return {};
    }
    return MySol;
  }
};

template <typename T>
bool TestEqualitySpannedSpaces(MyMatrix<T> const &M1, MyMatrix<T> const &M2) {
#ifdef SANITY_CHECK_MAT_MATRIX
  if (M1.cols() != M2.cols()) {
    std::cerr << "That case is actually not allowed\n";
    throw TerminalException{1};
  }
  if (RankMat(M1) != M1.rows()) {
    std::cerr << "M1 number of rows should match its rank\n";
    throw TerminalException{1};
  }
  if (RankMat(M2) != M2.rows()) {
    std::cerr << "M2 number of rows should match its rank\n";
    throw TerminalException{1};
  }
#endif
  // We make the assumption that the input is valid, that is
  // that the number of rows is equal to the rank.
  if (M1.rows() != M2.rows()) {
    return false;
  }
  SolutionMatRepetitive<T> smr(M1);
  for (int irow = 0; irow < M2.rows(); irow++) {
    MyVector<T> V2 = GetMatrixRow(M2, irow);
    std::optional<MyVector<T>> opt = smr.GetSolution(V2);
    if (!opt) {
      return false;
    }
  }
  return true;
}

/*
  We can actually do a little bit better for the solution to avoid repeating
  the preprocessing.
 */
template <typename T>
std::optional<MyMatrix<T>> ListSolutionMat(MyMatrix<T> const &eMat,
                                           MyMatrix<T> const &LVect) {
  int n_vect = LVect.rows();
  int dim = eMat.rows();
  MyMatrix<T> TheSol(n_vect, dim);
  for (int i_vect = 0; i_vect < n_vect; i_vect++) {
    MyVector<T> V = GetMatrixRow(LVect, i_vect);
    std::optional<MyVector<T>> opt = SolutionMat(eMat, V);
    if (!opt)
      return {};
    MyVector<T> const &V2 = *opt;
    AssignMatrixRow(TheSol, i_vect, V2);
  }
  return TheSol;
}

template <typename T>
MyMatrix<T> ExpressVectorsInIndependentFamilt(MyMatrix<T> const &VF,
                                              MyMatrix<T> const &IVF) {
  int n_vect = VF.rows();
  int dim = IVF.rows();
  MyMatrix<T> P(n_vect, dim);
  for (int i = 0; i < n_vect; i++) {
    MyVector<T> eV = GetMatrixRow(VF, i);
    std::optional<MyVector<T>> opt = SolutionMat(IVF, eV);
    if (!opt) {
      std::cerr << "VF : i=" << i << " not expressed in term of IVF\n";
      throw TerminalException{1};
    }
    AssignMatrixRow(P, i, *opt);
  }
  return P;
}

template <typename T>
MyVector<T> SolutionMat_LeastSquare(MyMatrix<T> const &M,
                                    MyVector<double> const &V) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  // Msqr should have size (nbCol, nbCol)
  MyMatrix<double> Msqr = M.transpose() * M;
  MyVector<double> B = ZeroVector<double>(nbCol);
  for (int j = 0; j < nbCol; j++)
    for (int i = 0; i < nbRow; i++)
      B(i, j) += V(i) * M(i, j);
  //
  // We have now a linear system to solve
  //
  Eigen::FullPivLU<MyMatrix<double>> solver;
  solver.compute(Msqr);
  MyVector<double> eSol = solver.solve(B);
  //
  return eSol;
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIXSOLUTIONMAT_H_
// clang-format on
