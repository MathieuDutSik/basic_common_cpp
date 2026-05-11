// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MATRIXTYPESSPARSE_H_
#define SRC_MATRIX_MATRIXTYPESSPARSE_H_

#include "MatrixTypes.h"
#include <Eigen/Sparse>

#if defined INCLUDE_NUMBER_THEORY_BOOST_GMP_INT ||                             \
    defined INCLUDE_NUMBER_THEORY_BOOST_CPP_INT
#include <boost/multiprecision/eigen.hpp>
#endif

#include <vector>

template <typename T>
using MySparseMatrix = Eigen::SparseMatrix<T, Eigen::ColMajor>;

namespace boost::serialization {

//
// MySparseMatrix data type
//

template <class Archive, typename T>
inline void load(Archive &ar, MySparseMatrix<T> &val,
                 [[maybe_unused]] const unsigned int version) {
  int nbRow, nbCol, nnz;
  ar &make_nvp("rows", nbRow);
  ar &make_nvp("cols", nbCol);
  ar &make_nvp("nnz", nnz);
  using T2 = Eigen::Triplet<T>;
  std::vector<T2> tripletList(nnz);
  for (int iNNZ = 0; iNNZ < nnz; iNNZ++) {
    int iRow, iCol;
    T eVal;
    ar &make_nvp("iRow", iRow);
    ar &make_nvp("iCol", iCol);
    ar &make_nvp("eVal", eVal);
    tripletList[iNNZ] = T2(iRow, iCol, eVal);
  }
  MySparseMatrix<T> SpMat(nbRow, nbCol);
  SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
  val = SpMat;
}

template <class Archive, typename T>
inline void save(Archive &ar, MySparseMatrix<T> const &val,
                 [[maybe_unused]] const unsigned int version) {
  int nbRow = val.rows();
  int nbCol = val.cols();
  int nnz = val.nonZeros();
  ar &make_nvp("rows", nbRow);
  ar &make_nvp("cols", nbCol);
  ar &make_nvp("nnz", nnz);
  for (int k = 0; k < val.outerSize(); ++k)
    for (typename MySparseMatrix<T>::InnerIterator it(val, k); it; ++it) {
      T eVal = it.value();
      // row index
      int iRow = it.row();
      // col index (here it is equal to k)
      int iCol = it.col();
      ar &make_nvp("iRow", iRow);
      ar &make_nvp("iCol", iCol);
      ar &make_nvp("eVal", eVal);
    }
}

template <class Archive, typename T>
inline void serialize(Archive &ar, MySparseMatrix<T> &val,
                      const unsigned int version) {
  split_free(ar, val, version);
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

template <typename T2, typename T1>
MySparseMatrix<T2> UniversalSparseMatrixConversion(MySparseMatrix<T1> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  int nnz = M.nonZeros();
  using Ttrip = Eigen::Triplet<T2>;
  std::vector<Ttrip> tripletList(nnz);
  int pos = 0;
  for (int k = 0; k < M.outerSize(); ++k)
    for (typename MySparseMatrix<T1>::InnerIterator it(M, k); it; ++it) {
      T1 eVal1 = it.value();
      T2 eVal2 = UniversalScalarConversion<T2,T1>(eVal1);
      // row index
      int iRow = it.row();
      // col index (here it is equal to k)
      int iCol = it.col();
      tripletList[pos] = Ttrip(iRow, iCol, eVal2);
      pos++;
    }
  MySparseMatrix<T2> SpMat = MySparseMatrix<T2>(nbRow, nbCol);
  SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
  return SpMat;
}

template <typename T>
MySparseMatrix<T> SparseMatrixSelectRows(MySparseMatrix<T> const &M, std::vector<int> const& l_rows) {
  int nbRow_orig = M.rows();
  int nbRow = l_rows.size();
  std::vector<int> map(nbRow_orig,-1);
  for (int i=0; i<nbRow; i++) {
    map[l_rows[i]] = i;
  }
  int nbCol = M.cols();
  using Ttrip = Eigen::Triplet<T>;
  std::vector<Ttrip> tripletList;
  for (int k = 0; k < M.outerSize(); ++k)
    for (typename MySparseMatrix<T>::InnerIterator it(M, k); it; ++it) {
      // row index
      int iRow_orig = it.row();
      int iRow = map[iRow_orig];
      if (iRow >= 0) {
        T eVal = it.value();
        int iCol = it.col();
        Ttrip ent(iRow, iCol, eVal);
        tripletList.push_back(ent);
      }
    }
  MySparseMatrix<T> SpMat = MySparseMatrix<T>(nbRow, nbCol);
  SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
  return SpMat;
}

template <typename T>
void WriteSparseMatrixGAP(std::ostream &os, MySparseMatrix<T> const &TheMat) {
  size_t nbRow = TheMat.rows();
  size_t nbCol = TheMat.cols();
  struct PairCV {
    int iCol;
    T eVal;
  };
  std::vector<std::vector<PairCV>> LLPair(nbRow);
  for (int k = 0; k < TheMat.outerSize(); ++k)
    for (typename MySparseMatrix<T>::InnerIterator it(TheMat, k); it; ++it) {
      T eVal = it.value();
      int iRow = it.row();
      int iCol = it.col();
      PairCV ePair{iCol, eVal};
      LLPair[iRow].push_back(ePair);
    }
  os << "rec(nbRow:=" << nbRow << ", nbCol:=" << nbCol << ", ListEntries:=[";
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    if (iRow > 0)
      os << ",\n";
    size_t len = LLPair[iRow].size();
    os << "rec(ListCol:=[";
    for (size_t i = 0; i < len; i++) {
      if (i > 0)
        os << ", ";
      int eCol = LLPair[iRow][i].iCol + 1;
      os << eCol;
    }
    os << "], ListVal:=[";
    for (size_t i = 0; i < len; i++) {
      if (i > 0)
        os << ", ";
      os << LLPair[iRow][i].eVal;
    }
    os << "])";
  }
  os << "])";
}

template <typename T> MySparseMatrix<T> ReadSparseMatrix(std::istream &is) {
  if (!is.good()) {
    std::cerr
        << "ReadSparseMatrix operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  T eVal;
  int nbRow, nbCol, nnz;
  int iRow, iCol;
  is >> nbRow >> nbCol >> nnz;
  using T2 = Eigen::Triplet<T>;
  std::vector<T2> tripletList(nnz);
  for (int iNNZ = 0; iNNZ < nnz; iNNZ++) {
    is >> iRow >> iCol >> eVal;
    tripletList[iNNZ] = T2(iRow, iCol, eVal);
    if (iRow >= nbRow) {
      std::cerr << "iRow is too large\n";
      std::cerr << " iRow=" << iRow << "\n";
      std::cerr << "nbRow=" << nbRow << "\n";
      throw TerminalException{1};
    }
    if (iCol >= nbCol) {
      std::cerr << "iCol is too large\n";
      std::cerr << " iCol=" << iCol << "\n";
      std::cerr << "nbCol=" << nbCol << "\n";
      throw TerminalException{1};
    }
  }
  MySparseMatrix<T> SpMat(nbRow, nbCol);
  SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
  return SpMat;
}

template <typename T>
void WriteSparseMatrix(std::ostream &os, MySparseMatrix<T> const &eMat) {
  int nbRow = eMat.rows();
  int nbCol = eMat.cols();
  int nnz = eMat.nonZeros();
  os << nbRow << " " << nbCol << " " << nnz;
  int nb = 0;
  for (int k = 0; k < eMat.outerSize(); ++k)
    for (typename MySparseMatrix<T>::InnerIterator it(eMat, k); it; ++it) {
      T eVal = it.value();
      int iRow = it.row();
      // col index (here it is equal to k)
      int iCol = it.col();
      os << iRow << " " << iCol << " " << eVal << "\n";
      nb++;
    }
  if (nb != nnz) {
    std::cerr << "Logical error in our understanding of \n";
    throw TerminalException{1};
  }
}

template <typename T>
MyMatrix<T> MyMatrixFromSparseMatrix(MySparseMatrix<T> const &eMat) {
  int nbRow = eMat.rows();
  int nbCol = eMat.cols();
  MyMatrix<T> M;
  M.setConstant(nbRow, nbCol, T(0));
  for (int k = 0; k < eMat.outerSize(); ++k)
    for (typename MySparseMatrix<T>::InnerIterator it(eMat, k); it; ++it) {
      T eVal = it.value();
      int iRow = it.row();
      // col index (here it is equal to k)
      int iCol = it.col();
      M(iRow, iCol) = eVal;
    }
  return M;
}

// clang-format off
#endif  // SRC_MATRIX_MATRIXTYPESSPARSE_H_
// clang-format on
