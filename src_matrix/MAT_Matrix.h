// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_MATRIX_H_
#define SRC_MATRIX_MAT_MATRIX_H_

// For reference over Eigen, see for example
// http://eigen.tuxfamily.org/dox/AsciiQuickReference.txt
//

// clang-format off
#include "Basic_file.h"
#include "Basic_string.h"
#include "MatrixTypes.h"
#include "Temp_common.h"
#include "Timings.h"
#include "hash_functions.h"
#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>
// clang-format on

template <typename T> struct is_mymatrix<MyMatrix<T>> {
  static const bool value = true;
};

namespace boost::serialization {

//
// MyVector data type
//

template <class Archive, typename T>
inline void serialize(Archive &ar, MyVector<T> &vect,
                      [[maybe_unused]] const unsigned int version) {
  int rows = vect.rows();
  ar &make_nvp("rows", rows);
  // matrix.resize is a no-op if size does not change!
  vect.resize(rows, 1);
  // always save/load row-major
  for (int r = 0; r < rows; ++r)
    ar &make_nvp("val", vect(r));
}

//
// MyMatrix data type
//

template <class Archive, typename T>
inline void serialize(Archive &ar, MyMatrix<T> &matrix,
                      [[maybe_unused]] const unsigned int version) {
  int rows = matrix.rows();
  int cols = matrix.cols();
  ar &make_nvp("rows", rows);
  ar &make_nvp("cols", cols);
  // matrix.resize is a no-op if size does not change!
  matrix.resize(rows, cols);
  // always save/load row-major
  for (int r = 0; r < rows; ++r)
    for (int c = 0; c < cols; ++c)
      ar &make_nvp("val", matrix(r, c));
}

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

template <typename T> bool IsIdentity(MyMatrix<T> const &M) {
  int len = M.rows();
  for (int i = 0; i < len; i++) {
    for (int j = 0; j < len; j++) {
      if (i == j && M(i, j) != 1)
        return false;
      if (i != j && M(i, j) != 0)
        return false;
    }
  }
  return true;
}

//
// Matrix type conversion
//

template <typename T2, typename T1>
MyVector<T2> UniversalVectorConversion(MyVector<T1> const &V) {
  int n = V.size();
  MyVector<T2> eRet(n);
  for (int i = 0; i < n; i++) {
    T1 eVal1 = V(i);
    T2 eVal2 = UniversalScalarConversion<T2, T1>(eVal1);
    eRet(i) = eVal2;
  }
  return eRet;
}

template <typename T2, typename T1>
std::optional<MyVector<T2>> UniversalVectorConversionCheck(MyVector<T1> const &V) {
  int n = V.size();
  MyVector<T2> V_ret(n);
  try {
    T2 ret;
    for (int i = 0; i < n; i++) {
      stc<T1> stc_a{V(i)};
      TYPE_CONVERSION(stc_a, ret);
      V_ret(i) = ret;
    }
  } catch (ConversionException &e) {
    return {};
  }
  return V_ret;
}

template <typename T2, typename T1>
MyMatrix<T2> UniversalMatrixConversion(MyMatrix<T1> const &M) {
  int n_rows = M.rows();
  int n_cols = M.cols();
  MyMatrix<T2> eRet(n_rows, n_cols);
  for (int j = 0; j < n_cols; j++)
    for (int i = 0; i < n_rows; i++)
      eRet(i, j) = UniversalScalarConversion<T2, T1>(M(i, j));
  return eRet;
}

template <typename T2, typename T1>
std::optional<MyMatrix<T2>> UniversalMatrixConversionCheck(MyMatrix<T1> const &M) {
  int n_rows = M.rows();
  int n_cols = M.cols();
  MyMatrix<T2> M_ret(n_rows, n_cols);
  try {
    T2 ret;
    for (int j = 0; j < n_cols; j++) {
      for (int i = 0; i < n_rows; i++) {
        stc<T1> stc_a{M(i,j)};
        TYPE_CONVERSION(stc_a, ret);
        M_ret(i, j) = ret;
      }
    }
  } catch (ConversionException &e) {
    return {};
  }
  return M_ret;
}

template <typename T2, typename T1>
std::vector<MyMatrix<T2>>
UniversalStdVectorMatrixConversion(std::vector<MyMatrix<T1>> const &ListM) {
  size_t n_mat = ListM.size();
  std::vector<MyMatrix<T2>> ListM_ret(n_mat);
  for (size_t i_mat = 0; i_mat < n_mat; i_mat++)
    ListM_ret[i_mat] = UniversalMatrixConversion<T2, T1>(ListM[i_mat]);
  return ListM_ret;
}

template <typename T>
MyVector<T> VectorWithIdenticalEntries(int const &len, T const &val) {
  MyVector<T> V(len);
  for (int i = 0; i < len; i++)
    V(i) = val;
  return V;
}

template <typename T1, typename T2>
MySparseMatrix<T2> ConvertSparseMatrix(MySparseMatrix<T1> const &M,
                                       std::function<T2(T1)> const &f) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  int nnz = M.nonZeros();
  using Ttrip = Eigen::Triplet<T2>;
  std::vector<Ttrip> tripletList(nnz);
  int nb = 0;
  for (int k = 0; k < M.outerSize(); ++k)
    for (typename MySparseMatrix<T1>::InnerIterator it(M, k); it; ++it) {
      T1 eVal1 = it.value();
      T2 eVal2 = f(eVal1);
      // row index
      int iRow = it.row();
      // col index (here it is equal to k)
      int iCol = it.col();
      tripletList[nb] = Ttrip(iRow, iCol, eVal2);
      nb++;
    }
  MySparseMatrix<T2> SpMat = MySparseMatrix<T2>(nbRow, nbCol);
  SpMat.setFromTriplets(tripletList.begin(), tripletList.end());
  return SpMat;
}

template <typename T>
MyMatrix<T> DropColumn(MyMatrix<T> const &M, int const &idx_drop) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  MyMatrix<T> Mred(nbRow, nbCol - 1);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    int pos = 0;
    for (int iCol = 0; iCol < nbCol; iCol++) {
      if (iCol != idx_drop) {
        Mred(iRow, pos) = M(iRow, iCol);
        pos++;
      }
    }
  }
  return Mred;
}

//
// ReadWriteMatrices and vectors
//

template <typename T> MyMatrix<T> ReadMatrix(std::istream &is) {
  if (!is.good()) {
    std::cerr << "ReadMatrix operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  T eVal;
  int nbRow, nbCol;
  is >> nbRow >> nbCol;
  if (nbRow < 0 || nbCol < 0) {
    std::cerr << "We should have nbRow > 0 and nbCol > 0\n";
    std::cerr << "But we have nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    throw TerminalException{1};
  }
  MyMatrix<T> TheMat(nbRow, nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++)
    for (int iCol = 0; iCol < nbCol; iCol++) {
      is >> eVal;
      TheMat(iRow, iCol) = eVal;
    }
  return TheMat;
}

template <typename T> MyMatrix<T> ReadMatrixFile(std::string const &file_name) {
  if (!IsExistingFile(file_name)) {
    std::cerr << "Error in ReadMatrixFile\n";
    std::cerr << "file_name=" << file_name << " does not appear to exist\n";
    throw TerminalException{1};
  }
  std::ifstream is(file_name);
  return ReadMatrix<T>(is);
}

template <typename T>
std::vector<MyMatrix<T>> ReadListMatrix(std::istream &is) {
  int n_mat;
  is >> n_mat;
  std::vector<MyMatrix<T>> ListMat;
  for (int i_mat = 0; i_mat < n_mat; i_mat++) {
    MyMatrix<T> eMat = ReadMatrix<T>(is);
    ListMat.push_back(eMat);
  }
  return ListMat;
}

template <typename T>
std::vector<MyMatrix<T>> ReadListMatrixFile(std::string const &eFile) {
  if (!IsExistingFile(eFile)) {
    std::cerr << "Error in ReadListMatrixFile\n";
    std::cerr << "eFile=" << eFile << " does not appear to exist\n";
    throw TerminalException{1};
  }
  std::ifstream is(eFile);
  return ReadListMatrix<T>(is);
}

template <typename T>
void WriteListMatrix(std::ostream &os,
                     std::vector<MyMatrix<T>> const &ListMat) {
  size_t n_mat = ListMat.size();
  os << n_mat << "\n";
  for (size_t i_mat = 0; i_mat < n_mat; i_mat++) {
    WriteMatrix(os, ListMat[i_mat]);
  }
}

template <typename T>
void WriteListMatrixFile(std::string const &eFile,
                         std::vector<MyMatrix<T>> const &ListMat) {
  std::ofstream os(eFile);
  WriteListMatrix(os, ListMat);
}

template <typename T> std::pair<bool, T> ReadMatrixInfo(std::istream &is) {
  T eVal;
  int nbRow, nbCol;
  is >> nbRow >> nbCol;
  if (nbRow < 0 || nbCol < 0) {
    std::cerr << "We should have nbRow > 0 and nbCol > 0\n";
    std::cerr << "But we have nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    throw TerminalException{1};
  }
  bool is_integral = true;
  T TheMax = 0;
  for (int iRow = 0; iRow < nbRow; iRow++)
    for (int iCol = 0; iCol < nbCol; iCol++) {
      is >> eVal;
      if (!IsInteger(eVal))
        is_integral = false;
      T eValAbs = T_abs(eVal);
      if (eValAbs > TheMax)
        TheMax = eValAbs;
    }
  return {is_integral, TheMax};
}

template <typename T> MyMatrix<T> ReadMatrixLrsCdd(std::istream &is) {
  if (!is.good()) {
    std::cerr << "ReadMatrixLrs operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  std::string header;
  is >> header;
  if (header != "V-representation" && header != "H-representation") {
    std::cerr << "Error while reading header\n";
    std::cerr << "header = " << header << "\n";
    throw TerminalException{1};
  }
  //
  std::string head2;
  is >> head2;
  if (head2 != "begin") {
    std::cerr << "Error while reading\n";
    std::cerr << "head2 = " << head2 << "\n";
    throw TerminalException{1};
  }
  //
  int nbRow, nbCol;
  std::string str;
  is >> nbRow >> nbCol >> str;
  if (str != "integer" && str != "rational" && str != "double") {
    std::cerr << "Error while reading\n";
    std::cerr << "str = " << str << "\n";
    std::cerr << "But it should be integer, rational or double. Note that it "
                 "is not being used\n";
    throw TerminalException{1};
  }
  MyMatrix<T> TheMat(nbRow, nbCol);
  for (int iRow = 0; iRow < nbRow; iRow++)
    for (int iCol = 0; iCol < nbCol; iCol++) {
      T eVal;
      is >> eVal;
      TheMat(iRow, iCol) = eVal;
    }
  //
  std::string footer;
  is >> footer;
  if (footer != "end") {
    std::cerr << "Error while reading\n";
    std::cerr << "footer = " << footer << " while it should be \"end\"\n";
    std::cerr << "What has been read (if that helps you):\n";
    std::cerr << "nbRow=" << nbRow << " nbCol=" << nbCol << "\n";
    for (int iRow = 0; iRow < nbRow; iRow++) {
      std::cerr << "iRow=" << iRow << " :";
      for (int iCol = 0; iCol < nbCol; iCol++)
        std::cerr << " " << TheMat(iRow, iCol);
      std::cerr << "\n";
    }
    throw TerminalException{1};
  }
  return TheMat;
}

template <typename T1, typename T2>
bool IsEqualSizeMatrices(MyMatrix<T1> const &X, MyMatrix<T2> const &Y) {
  if (X.rows() != Y.rows())
    return false;
  if (X.cols() != Y.cols())
    return false;
  return true;
}

template <typename T> std::string StringSizeMatrix(MyMatrix<T> const &X) {
  int nbRow = X.rows();
  int nbCol = X.cols();
  return std::to_string(nbRow) + " / " + std::to_string(nbCol);
}

template <typename T> std::string MinMaxMatrix(MyMatrix<T> const &X) {
  T minV = X.minCoeff();
  T maxV = X.maxCoeff();
  return std::string("min/max=") + std::to_string(minV) + " / " +
         std::to_string(maxV);
}

template <typename T> MyVector<T> ReadVector(std::istream &is) {
  if (!is.good()) {
    std::cerr << "ReadVector operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  T eVal;
  size_t nbRow;
  is >> nbRow;
  MyVector<T> eVect = MyVector<T>(nbRow);
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    is >> eVal;
    eVect(iRow) = eVal;
  }
  return eVect;
}

template <typename T> MyVector<T> ReadVectorFile(std::string const &file_name) {
  if (!IsExistingFile(file_name)) {
    std::cerr << "Error in ReadVectorFile\n";
    std::cerr << "file_name=" << file_name << " does not appear to exist\n";
    throw TerminalException{1};
  }
  std::ifstream is(file_name);
  return ReadVector<T>(is);
}

template <typename T> std::vector<T> ReadStdVector(std::istream &is) {
  if (!is.good()) {
    std::cerr << "ReadStdVector operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  T eVal;
  size_t nbRow;
  is >> nbRow;
  std::vector<T> eVect(nbRow);
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    is >> eVal;
    eVect[iRow] = eVal;
  }
  return eVect;
}

template <typename T>
void WriteMatrix(std::ostream &os, MyMatrix<T> const &TheMat) {
  size_t nbRow = TheMat.rows();
  size_t nbCol = TheMat.cols();
  //  TerminalEnding();
  os << nbRow << " " << nbCol << "\n";
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      os << " " << TheMat(iRow, iCol);
    os << "\n";
  }
}

template <typename T>
void WriteMatrixFile(std::string const &eFile, MyMatrix<T> const &TheMat) {
  std::ofstream os(eFile);
  WriteMatrix(os, TheMat);
}

template <typename T>
void WriteMatrixNice(std::ostream &os, MyMatrix<T> const &M) {
  size_t nbRow = M.rows();
  size_t nbCol = M.cols();
  //  TerminalEnding();
  os << nbRow << " " << nbCol << "\n";
  std::vector<std::vector<std::string>> LLStr;
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    std::vector<std::string> LStr;
    for (size_t iCol = 0; iCol < nbCol; iCol++) {
      std::stringstream s;
      s << M(iRow, iCol);
      std::string converted(s.str());
      LStr.push_back(converted);
    }
    LLStr.push_back(LStr);
  }
  std::vector<size_t> l_max_nchar(nbCol);
  for (size_t iCol = 0; iCol < nbCol; iCol++) {
    size_t max_nchar = 0;
    for (size_t iRow = 0; iRow < nbRow; iRow++)
      max_nchar = std::max(max_nchar, LLStr[iRow][iCol].size());
    l_max_nchar[iCol] = max_nchar;
  }
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    for (size_t iCol = 0; iCol < nbCol; iCol++) {
      std::string str = LLStr[iRow][iCol];
      size_t n_sp = l_max_nchar[iCol] - str.size();
      os << " " << str;
      for (size_t i = 0; i < n_sp; i++)
        os << " ";
    }
    os << "\n";
  }
}

template <typename T>
void WriteMatrixMatlab(std::ostream &os, MyMatrix<T> const &TheMat) {
  int nbRow = TheMat.rows();
  int nbCol = TheMat.cols();
  for (int iRow = 0; iRow < nbRow; iRow++) {
    for (int iCol = 0; iCol < nbCol; iCol++)
      os << " " << TheMat(iRow, iCol);
    os << "\n";
  }
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

template <typename T>
void WriteMatrixPYTHON(std::ostream &os, MyMatrix<T> const &TheMat) {
  size_t nbRow = TheMat.rows();
  size_t nbCol = TheMat.cols();
  os << "[";
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    if (iRow > 0) {
      os << ",";
    }
    os << "[";
    for (size_t iCol = 0; iCol < nbCol; iCol++) {
      T eVal = TheMat(iRow, iCol);
      if (iCol > 0)
        os << ",";
      os << eVal;
    }
    os << "]";
  }
  os << "]";
}

template <typename T>
std::string StringMatrixPYTHON(MyMatrix<T> const& TheMat) {
  std::ostringstream os;
  WriteMatrixPYTHON(os, TheMat);
  return os.str();
}

template <typename T>
void WriteMatrixGAP_gen(std::ostream &os, MyMatrix<T> const &TheMat,
                        bool const &as_line) {
  size_t nbRow = TheMat.rows();
  size_t nbCol = TheMat.cols();
  os << "[ ";
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    if (iRow > 0) {
      os << ",";
      if (!as_line)
        os << "\n";
    }
    os << "[";
    for (size_t iCol = 0; iCol < nbCol; iCol++) {
      T eVal = TheMat(iRow, iCol);
      if (iCol > 0)
        os << ",";
      os << " " << eVal;
    }
    os << " ]";
  }
  os << " ]";
}

template <typename T>
void WriteMatrixGAP(std::ostream &os, MyMatrix<T> const &TheMat) {
  WriteMatrixGAP_gen(os, TheMat, false);
}

template <typename T> std::string StringMatrixGAP(MyMatrix<T> const &TheMat) {
  std::ostringstream os;
  WriteMatrixGAP(os, TheMat);
  return os.str();
}

template <typename T>
void WriteListMatrixGAP(std::ostream &os,
                        std::vector<MyMatrix<T>> const &ListMat) {
  os << "[";
  bool IsFirst = true;
  for (auto &eMat : ListMat) {
    if (!IsFirst)
      os << ",\n";
    IsFirst = false;
    WriteMatrixGAP(os, eMat);
  }
  os << "]";
}

template<typename T>
void WriteListMatrixPYTHON(std::ostream& os,
                           std::vector<MyMatrix<T>> const &ListMat) {
  os << "[";
  bool IsFirst = true;
  for (auto &eMat : ListMat) {
    if (!IsFirst)
      os << ",";
    IsFirst = false;
    WriteMatrixPYTHON(os, eMat);
  }
  os << "]";
}

template<typename T>
std::string StringListMatrixPYTHON(std::vector<MyMatrix<T>> const &ListMat) {
  std::ostringstream os;
  WriteListMatrixPYTHON(os, ListMat);
  return os.str();
}

template <typename T>
void WriteListMatrixFileGAP(std::string &eFile,
                            std::vector<MyMatrix<T>> const &ListMat) {
  std::ofstream os_out(eFile);
  os_out << "return ";
  WriteListMatrixGAP(os_out, ListMat);
  os_out << ";\n";
}

template <typename T>
std::string StringMatrixGAP_line(MyMatrix<T> const &TheMat) {
  std::ostringstream os;
  WriteMatrixGAP_gen(os, TheMat, true);
  return os.str();
}

template <typename T>
void WriteVectorMatrixGAP(std::ostream &os,
                          std::vector<MyMatrix<T>> const &l_mat) {
  os << "[";
  for (size_t i = 0; i < l_mat.size(); i++) {
    if (i > 0)
      os << ",";
    WriteMatrixGAP(os, l_mat[i]);
  }
  os << "]";
}

template <typename T>
void WriteVectorMatrixPYTHON(std::ostream &os,
                             std::vector<MyMatrix<T>> const &l_mat) {
  os << "[";
  for (size_t i = 0; i < l_mat.size(); i++) {
    if (i > 0)
      os << ",";
    WriteMatrixPYTHON(os, l_mat[i]);
  }
  os << "]";
}

template <typename T>
void WriteMatrixGAPfile(std::string const &eFile, MyMatrix<T> const &TheMat) {
  std::ofstream os(eFile);
  os << "return ";
  WriteMatrixGAP(os, TheMat);
  os << ";\n";
}

template <typename T>
std::ostream &operator<<(std::ostream &os, MyVector<T> const &TheVec) {
  int n = TheVec.size();
  if (n > 0) {
    os << TheVec(0);
    for (int i = 1; i < n; i++)
      os << " " << TheVec(i);
  }
  return os;
}

template <typename T>
void WriteVector(std::ostream &os, MyVector<T> const &TheVec) {
  int n = TheVec.size();
  os << n << "\n";
  for (int i = 0; i < n; i++)
    os << " " << TheVec(i);
  os << "\n";
}

template <typename T>
void WriteVectorNoDim(std::ostream &os, MyVector<T> const &TheVec) {
  int n = TheVec.size();
  for (int i = 0; i < n; i++)
    os << " " << TheVec(i);
  os << "\n";
}

template <typename T>
void WriteVectorPYTHON(std::ostream &os, MyVector<T> const &TheVec) {
  int n = TheVec.size();
  os << "[";
  for (int i = 0; i < n; i++) {
    if (i > 0)
      os << ",";
    os << TheVec(i);
  }
  os << "]";
}

template <typename T>
std::string StringVectorPYTHON(MyVector<T> const &TheVec) {
  std::ostringstream os;
  WriteVectorPYTHON(os, TheVec);
  return os.str();
}

template <typename T>
void WriteVectorGAP(std::ostream &os, MyVector<T> const &TheVec) {
  int n = TheVec.size();
  os << "[ ";
  for (int i = 0; i < n; i++) {
    if (i > 0)
      os << ", ";
    os << TheVec(i);
  }
  os << " ]";
}

template <typename T>
void WriteVectorFile(std::string const &eFile, MyVector<T> const &TheV) {
  std::ofstream os(eFile);
  WriteVector(os, TheV);
}

template <typename T> std::string StringVector(MyVector<T> const &TheVec) {
  std::ostringstream os;
  int n = TheVec.size();
  for (int i = 0; i < n; i++)
    os << " " << TheVec(i);
  return os.str();
}

template <typename T> std::string StringVectorGAP(MyVector<T> const &TheVec) {
  std::ostringstream os;
  WriteVectorGAP(os, TheVec);
  return os.str();
}

template <typename T> T VectorSum(MyVector<T> const &eVect) {
  T eSum(0);
  int siz = eVect.size();
  for (int i = 0; i < siz; i++)
    eSum += eVect(i);
  return eSum;
}

template <typename T>
T ScalarProduct(MyVector<T> const &V1, MyVector<T> const &V2) {
  if (V1.size() != V2.size()) {
    std::cerr << "Vectors of wrong sizes\n";
    throw TerminalException{1};
  }
  size_t siz = V1.size();
  T eSum(0);
  for (size_t i = 0; i < siz; i++)
    eSum += V1(i) * V2(i);
  return eSum;
}

// compute the scalar product of V with the lines of eMat
template <typename T>
MyVector<T> ListScalarProduct(MyVector<T> const &V, MyMatrix<T> const &eMat) {
  int dim = V.size();
  int nbRow = eMat.rows();
  int nbCol = eMat.cols();
  if (nbCol != dim) {
    std::cerr << "dim should equal nbCol\n";
    throw TerminalException{1};
  }
  MyVector<T> retVect(nbRow);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    T eSum(0);
    for (int iCol = 0; iCol < dim; iCol++)
      eSum += eMat(iRow, iCol) * V(iCol);
    retVect(iRow) = eSum;
  }
  return retVect;
}

template <typename T>
MyMatrix<T> ZeroMatrix(int const &nbRow, int const &nbCol) {
  MyMatrix<T> retMat;
  retMat.setConstant(nbRow, nbCol, T(0));
  return retMat;
}

template <typename T>
MyMatrix<T> ConstantMatrix(int const &nbRow, int const &nbCol, T const &eVal) {
  MyMatrix<T> retMat;
  retMat.setConstant(nbRow, nbCol, eVal);
  return retMat;
}

template <typename T>
MyMatrix<T> TranspositionMatrix(int const &n, int const &i, int const &j) {
  MyMatrix<T> M;
  T val(0);
  M.setConstant(n, n, val);
  for (int iCol = 0; iCol < n; iCol++) {
    int iRow = iCol;
    if (iCol == i)
      iRow = j;
    if (iCol == j)
      iRow = i;
    M(iRow, iCol) = 1;
  }
  return M;
}

template <typename T> MyVector<T> ZeroVector(int const &nbRow) {
  MyVector<T> retVect(nbRow);
  T eZero(0);
  for (int iRow = 0; iRow < nbRow; iRow++)
    retVect(iRow) = eZero;
  return retVect;
}

template <typename T> void TVec_ZeroAssignation(MyVector<T> &TheVect) {
  T eZero(0);
  for (int i = 0; i < TheVect.size(); i++)
    TheVect(i) = eZero;
}

template <typename T> void ZeroAssignation(MyMatrix<T> &TheMat) {
  int nbRow = TheMat.rows();
  int nbCol = TheMat.cols();
  T eVal(0);
  for (int iRow = 0; iRow < nbRow; iRow++)
    for (int iCol = 0; iCol < nbCol; iCol++)
      TheMat(iRow, iCol) = eVal;
}

template <typename T> MyMatrix<T> TransposedMat(MyMatrix<T> const &TheMat) {
  int nbCol = TheMat.cols();
  int nbRow = TheMat.rows();
  MyMatrix<T> TheTrans(nbCol, nbRow);
  for (int iCol = 0; iCol < nbCol; iCol++)
    for (int iRow = 0; iRow < nbRow; iRow++)
      TheTrans(iCol, iRow) = TheMat(iRow, iCol);
  return TheTrans;
}

// Compute the product XM
template <typename T>
MyVector<T> ProductVectorMatrix(MyVector<T> const &X, MyMatrix<T> const &M) {
  int nbCol = M.cols();
  int nbRow = M.rows();
  if (X.size() != nbRow) {
    std::cerr << "Error in the product X A\n";
    throw TerminalException{1};
  }
  MyVector<T> Vret(nbCol);
  for (int iCol = 0; iCol < nbCol; iCol++) {
    T sum(0);
    for (int iRow = 0; iRow < nbRow; iRow++)
      sum += M(iRow, iCol) * X(iRow);
    Vret(iCol) = sum;
  }
  return Vret;
}

template <typename T> MyMatrix<T> RankOneMatrix(MyVector<T> const &V) {
  int n = V.size();
  MyMatrix<T> retMat(n, n);
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++)
      retMat(i, j) = V(i) * V(j);
  return retMat;
}

// Use the companion matrix for computing the roots of a
// polynomial.
template <typename T>
std::vector<T> RootPolynomial(std::vector<T> const &eVect) {
  int nbEnt = eVect.size();
  MyMatrix<T> CompMat = ZeroMatrix<T>(nbEnt, nbEnt);
  for (int i = 0; i < nbEnt - 1; i++)
    CompMat(i + 1, i) = 1;
  for (int i = 0; i < nbEnt; i++)
    CompMat(i, nbEnt - 1) = -eVect[i];
  Eigen::EigenSolver<MyMatrix<T>> eig(CompMat);
  MyVector<T> ListEig = eig.eigenvalues();
  std::vector<T> ListRoot(nbEnt);
  for (int i = 0; i < nbEnt; i++)
    ListRoot[i] = ListEig(i);
  return ListRoot;
}

// form the product x*A with x a line vector and A a matrix
template <typename T>
MyVector<T> VectorMatrix(MyVector<T> const &eVect, MyMatrix<T> const &eMat) {
  int nbCol = eMat.cols();
  int nbRow = eMat.rows();
#ifdef SANITY_CHECK_MAT_MATRIX
  int n = eVect.size();
  if (n != nbRow) {
    std::cerr << "n should be equal to nbRow\n";
    throw TerminalException{1};
  }
#endif
  MyVector<T> rVect(nbCol);
  for (int iCol = 0; iCol < nbCol; iCol++) {
    T eSum = 0;
    for (int iRow = 0; iRow < nbRow; iRow++)
      eSum += eMat(iRow, iCol) * eVect(iRow);
    rVect(iCol) = eSum;
  }
  return rVect;
}

template <typename T>
void AssignMatrixRow(MyMatrix<T> &eMat, int const &iRow,
                     MyVector<T> const &eVect) {
  int nbCol = eMat.cols();
#ifdef SANITY_CHECK_MAT_MATRIX
  int n = eVect.size();
  if (n != nbCol) {
    std::cerr << "We should have eVect.size() = eMat.cols()\n";
    std::cerr << "eVect.size()=" << n << " eMat.cols()=" << nbCol << "\n";
    throw TerminalException{1};
  }
  int nbRow = eMat.rows();
  if (iRow >= nbRow) {
    std::cerr << "We should have iRow < nbRow\n";
    std::cerr << "iRow=" << iRow << " nbRow=" << nbRow << "\n";
    throw TerminalException{1};
  }
#endif
  for (int iCol = 0; iCol < nbCol; iCol++)
    eMat(iRow, iCol) = eVect(iCol);
}

template <typename T>
void AssignMatrixCol(MyMatrix<T> &eMat, int const &iCol,
                     MyVector<T> const &eVect) {
  int nbRow = eMat.rows();
  for (int iRow = 0; iRow < nbRow; iRow++)
    eMat(iRow, iCol) = eVect(iRow);
}

template <typename T>
bool IsVectorPositiveMultiple(MyVector<T> const &eVec1,
                              MyVector<T> const &eVec2) {
  int n = eVec1.size();
  bool IsAssign = false;
  T val2save, val1save;
  for (int i = 0; i < n; i++) {
    T eVal1 = eVec1(i);
    T eVal2 = eVec2(i);
    if (!IsAssign) {
      if (eVal1 != 0) {
        if (eVal2 == 0)
          return false;
        val2save = eVal2;
        val1save = eVal1;
        IsAssign = true;
      } else {
        if (eVal2 != 0)
          return false;
      }
    } else {
      T eDiff = val2save * eVal1 - val1save * eVal2;
      if (eDiff != 0)
        return false;
    }
  }
  if (!IsAssign) {
    std::cerr << "Vectors eVec1 is 0!\n";
    throw TerminalException{1};
  }
  if (val2save > 0 && val1save > 0)
    return true;
  if (val2save < 0 && val1save < 0)
    return true;
  return false;
}

template <typename T>
bool IsVectorMultiple(MyVector<T> const &eVec1, MyVector<T> const &eVec2) {
  int n = eVec1.size();
  bool IsAssign = false;
  T val2save, val1save;
  for (int i = 0; i < n; i++) {
    T eVal1 = eVec1(i);
    T eVal2 = eVec2(i);
    if (!IsAssign) {
      if (eVal1 != 0) {
        if (eVal2 == 0)
          return false;
        val2save = eVal2;
        val1save = eVal1;
        IsAssign = true;
      } else {
        if (eVal2 != 0)
          return false;
      }
    } else {
      T eDiff = val2save * eVal1 - val1save * eVal2;
      if (eDiff != 0)
        return false;
    }
  }
  if (!IsAssign) {
    std::cerr << "Vectors eVec1 is 0!\n";
    throw TerminalException{1};
  }
  return true;
}

template <typename T>
T DivideVector(MyVector<T> const &V1, MyVector<T> const &V2) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in DivideVector");
  int n = V1.size();
  for (int i = 0; i < n; i++)
    if (V1(i) != 0)
      return V1(i) / V2(i);
  std::cerr << "Error in DivideVector\n";
  throw TerminalException{1};
}

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
    iColFound = -1;
    prov1 = 0;
    for (iCol = iRow; iCol < nbCol; iCol++)
      if (iColFound == -1) {
        prov1 = Input(iRow, iCol);
        if (prov1 != 0) {
          iColFound = iCol;
          prov1 = 1 / prov1;
        }
      }
#ifdef SANITY_CHECK_MAT_MATRIX
    if (prov1 == 0) {
      std::cerr << "Error during the computation of the matrix inverse\n";
      throw TerminalException{1};
    }
#endif
    for (iRowB = 0; iRowB < nbRow; iRowB++)
      Output(iRowB, iColFound) *= prov1;
    for (iRowB = iRow; iRowB < nbRow; iRowB++)
      Input(iRowB, iColFound) *= prov1;
    for (iCol = 0; iCol < nbCol; iCol++)
      if (iCol != iColFound) {
        prov1 = Input(iRow, iCol);
        if (prov1 != 0) {
          for (iRowB = 0; iRowB < nbRow; iRowB++)
            Output(iRowB, iCol) -= prov1 * Output(iRowB, iColFound);
          for (iRowB = iRow; iRowB < nbRow; iRowB++)
            Input(iRowB, iCol) -= prov1 * Input(iRowB, iColFound);
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

template <typename T>
inline typename std::enable_if<is_ring_field<T>::value, MyMatrix<T>>::type
Inverse(MyMatrix<T> const &Input) {
  return InverseKernel(Input);
}

template <typename T>
inline typename std::enable_if<!is_ring_field<T>::value, MyMatrix<T>>::type
Inverse(MyMatrix<T> const &Input) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> InputF = UniversalMatrixConversion<Tfield, T>(Input);
  MyMatrix<Tfield> OutputF = InverseKernel(InputF);
  return UniversalMatrixConversion<T, Tfield>(OutputF);
}

/* This function is for rank calculation.
   Of course, it can be used for many other purpose:
   1> Selecting specific sets of rows and columns for reduction
   *   2> Computing a set of generators of the zero set.
   Initial matrix is of the form
   Input (nbRow, nbCol)
   The zero matrix is of the form
   NSP (dimKer, nbCol)
*/
template <typename T> struct SelectionRowCol {
  size_t TheRank;
  MyMatrix<T> NSP;
  std::vector<int> ListColSelect;
  std::vector<int> ListRowSelect;
};

template <typename T> struct RankTool {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in RankTool");

  RankTool(int const &eDim) : rank(0), dim(eDim), ListICol({}), ListVect({}) {}

  int insertion_operation(MyVector<T> &V) {
    for (int i_line = 0; i_line < rank; i_line++) {
      int eCol = ListICol[i_line];
      V -= V[eCol] * ListVect[i_line];
    }
    for (int i_dim = 0; i_dim < dim; i_dim++)
      if (V[i_dim] != 0) {
        V /= V[i_dim];
        return i_dim;
      }
    return -1;
  }

  void insert_if_indep(MyVector<T> &V) {
    int i_col = insertion_operation(V);
    if (i_col != -1) {
      rank++;
      ListICol.push_back(i_col);
      ListVect.push_back(V);
    }
  }

  bool is_in_span(MyVector<T> &V) {
    int i_col = insertion_operation(V);
    return i_col == -1;
  }

  int get_rank() const { return rank; }

private:
  int rank;
  int dim;
  std::vector<int> ListICol;
  std::vector<MyVector<T>> ListVect;
};

// The NSP array is assigned to NullspaceMat(TransposedMat(Input))
// The function requires the matrix type to be a field.
// After reflection this way is unavoidable as working only
// with ring operation would complicate the code quite a lot
// and likely lead to explosion of coefficient.
template <typename T, typename F>
SelectionRowCol<T> TMat_SelectRowCol_Kernel(size_t nbRow, size_t nbCol, F f) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in TMat_SelectRowCol_Kernel");
  size_t maxRank = nbRow;
  if (nbCol < maxRank)
    maxRank = nbCol;
  size_t sizMat = maxRank + 1;
  MyMatrix<T> provMat(sizMat, nbCol);
  std::vector<int> ListColSelect;
  std::vector<int> ListRowSelect;
  std::vector<uint8_t> ListColSelect01(nbCol, 0);
  size_t eRank = 0;
  size_t miss_val = std::numeric_limits<size_t>::max();
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    f(provMat, eRank, iRow);
    for (size_t iRank = 0; iRank < eRank; iRank++) {
      int eCol = ListColSelect[iRank];
      T eVal1 = provMat(eRank, eCol);
      if (eVal1 != 0) {
        for (size_t iCol = 0; iCol < nbCol; iCol++)
          provMat(eRank, iCol) -= eVal1 * provMat(iRank, iCol);
      }
    }
    auto get_firstnonzerocol_iife = [&]() -> size_t {
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        if (provMat(eRank, iCol) != 0)
          return iCol;
      return miss_val;
    };
    size_t FirstNonZeroCol = get_firstnonzerocol_iife();
    if (FirstNonZeroCol != miss_val) {
      ListColSelect.push_back(FirstNonZeroCol);
      ListRowSelect.push_back(iRow);
      ListColSelect01[size_t(FirstNonZeroCol)] = 1;
      T eVal2 = 1 / provMat(eRank, FirstNonZeroCol);
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        provMat(eRank, iCol) *= eVal2;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        T eVal1 = provMat(iRank, FirstNonZeroCol);
        if (eVal1 != 0) {
          for (size_t iCol = 0; iCol < nbCol; iCol++)
            provMat(iRank, iCol) -= eVal1 * provMat(eRank, iCol);
        }
      }
      eRank++;
    }
  }
  size_t nbVectZero = nbCol - eRank;
  MyMatrix<T> NSP = ZeroMatrix<T>(nbVectZero, nbCol);
  size_t nbVect = 0;
  for (size_t iCol = 0; iCol < nbCol; iCol++)
    if (ListColSelect01[iCol] == 0) {
      NSP(nbVect, iCol) = 1;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        int eCol = ListColSelect[iRank];
        NSP(nbVect, eCol) = -provMat(iRank, iCol);
      }
      nbVect++;
    }
  return {eRank, std::move(NSP), std::move(ListColSelect),
          std::move(ListRowSelect)};
}

// The NSP array is assigned to NullspaceMat(TransposedMat(Input))
// The function requires the matrix type to be a field.
// This version uses the maximum pivot for selecting the pivot
// which indicates that it should be used for floating point
// arithmetics. The pivot has to be aboved a selected threshold
// in order to be considered valid.
template <typename T, typename F>
SelectionRowCol<T> TMat_SelectRowColMaxPivot_Kernel(size_t nbRow, size_t nbCol,
                                                    F f, T const &threshold) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in TMat_SelectRowColMaxPivot_Kernel");
  size_t maxRank = nbRow;
  if (nbCol < maxRank)
    maxRank = nbCol;
  size_t sizMat = maxRank + 1;
  MyMatrix<T> provMat(sizMat, nbCol);
  std::vector<int> ListColSelect;
  std::vector<int> ListRowSelect;
  std::vector<uint8_t> ListColSelect01(nbCol, 0);
  size_t eRank = 0;
  size_t miss_val = std::numeric_limits<size_t>::max();
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    f(provMat, eRank, iRow);
    for (size_t iRank = 0; iRank < eRank; iRank++) {
      int eCol = ListColSelect[iRank];
      T eVal1 = provMat(eRank, eCol);
      if (eVal1 != 0) {
        for (size_t iCol = 0; iCol < nbCol; iCol++)
          provMat(eRank, iCol) -= eVal1 * provMat(iRank, iCol);
      }
    }
    auto get_firstnonzerocol_iife = [&]() -> size_t {
      T MaxPivot = 0;
      size_t selected_pivot = miss_val;
      for (size_t iCol = 0; iCol < nbCol; iCol++) {
        T abs_pivot = T_abs(provMat(eRank, iCol));
        if (abs_pivot > threshold && abs_pivot > MaxPivot) {
          MaxPivot = abs_pivot;
          selected_pivot = iCol;
        }
      }
      return selected_pivot;
    };
    size_t FirstNonZeroCol = get_firstnonzerocol_iife();
    if (FirstNonZeroCol != miss_val) {
      ListColSelect.push_back(FirstNonZeroCol);
      ListRowSelect.push_back(iRow);
      ListColSelect01[size_t(FirstNonZeroCol)] = 1;
      T eVal2 = 1 / provMat(eRank, FirstNonZeroCol);
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        provMat(eRank, iCol) *= eVal2;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        T eVal1 = provMat(iRank, FirstNonZeroCol);
        if (eVal1 != 0) {
          for (size_t iCol = 0; iCol < nbCol; iCol++)
            provMat(iRank, iCol) -= eVal1 * provMat(eRank, iCol);
        }
      }
      eRank++;
    }
  }
  size_t nbVectZero = nbCol - eRank;
  MyMatrix<T> NSP = ZeroMatrix<T>(nbVectZero, nbCol);
  size_t nbVect = 0;
  for (size_t iCol = 0; iCol < nbCol; iCol++)
    if (ListColSelect01[iCol] == 0) {
      NSP(nbVect, iCol) = 1;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        int eCol = ListColSelect[iRank];
        NSP(nbVect, eCol) = -provMat(iRank, iCol);
      }
      nbVect++;
    }
  return {eRank, std::move(NSP), std::move(ListColSelect),
          std::move(ListRowSelect)};
}

template <typename T>
MyMatrix<T> NullspaceMatSingleVector(MyVector<T> const &V) {
  int n = V.size();
  for (int i = 0; i < n; i++) {
    if (V(i) != 0) {
      MyMatrix<T> NSP = ZeroMatrix<T>(n - 1, n);
      int pos = 0;
      for (int j = 0; j < n; j++) {
        if (i != j) {
          NSP(pos, i) = -V(j);
          NSP(pos, j) = V(i);
          pos++;
        }
      }
      return NSP;
    }
  }
  std::cerr << "Failed to find a non-zero index\n";
  throw TerminalException{1};
}

template <typename T>
SelectionRowCol<T> TMat_SelectRowColMaxPivot(MyMatrix<T> const &Input,
                                             T const &threshold) {
  size_t nbRow = Input.rows();
  size_t nbCol = Input.cols();
  auto f = [&](MyMatrix<T> &M, size_t eRank, size_t iRow) -> void {
    M.row(eRank) = Input.row(iRow);
  };
  return TMat_SelectRowCol_Kernel<T>(nbRow, nbCol, f, threshold);
}

template <typename T>
SelectionRowCol<T>
TMat_SelectRowColMaxPivot_subset(MyMatrix<T> const &Input,
                                 std::vector<size_t> const &Vsubset,
                                 T const &threshold) {
  size_t nbRow = Vsubset.size();
  size_t nbCol = Input.cols();
  auto f = [&](MyMatrix<T> &M, size_t eRank, size_t iRow) -> void {
    M.row(eRank) = Input.row(Vsubset[iRow]);
  };
  return TMat_SelectRowCol_Kernel<T>(nbRow, nbCol, f, threshold);
}

template <typename T>
bool IsVectorInSpace(const SelectionRowCol<T> &eSelect, const MyVector<T> &V) {
  size_t dim_nsp = eSelect.NSP.rows();
  size_t n_cols = eSelect.NSP.cols();
  for (size_t i_nsp = 0; i_nsp < dim_nsp; i_nsp++) {
    T scal(0);
    for (size_t i_col = 0; i_col < n_cols; i_col++)
      scal += eSelect.NSP(i_nsp, i_col) * V(i_col);
    if (scal != 0)
      return false;
  }
  return true;
}

template <typename T>
SelectionRowCol<T> TMat_SelectRowCol(MyMatrix<T> const &Input) {
  size_t nbRow = Input.rows();
  size_t nbCol = Input.cols();
  auto f = [&](MyMatrix<T> &M, size_t eRank, size_t iRow) -> void {
    M.row(eRank) = Input.row(iRow);
  };
  return TMat_SelectRowCol_Kernel<T>(nbRow, nbCol, f);
}

template <typename T>
SelectionRowCol<T>
TMat_SelectRowCol_subset(MyMatrix<T> const &Input,
                         std::vector<size_t> const &Vsubset) {
  size_t nbRow = Vsubset.size();
  size_t nbCol = Input.cols();
  auto f = [&](MyMatrix<T> &M, size_t eRank, size_t iRow) -> void {
    M.row(eRank) = Input.row(Vsubset[iRow]);
  };
  return TMat_SelectRowCol_Kernel<T>(nbRow, nbCol, f);
}

template <typename T, typename F>
MyMatrix<T> NullspaceTrMat_Kernel(size_t nbRow, size_t nbCol, F f) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in NullspaceTrMat_Kernel");
  size_t maxRank = nbRow;
  if (nbCol < maxRank)
    maxRank = nbCol;
  size_t sizMat = maxRank + 1;
  MyMatrix<T> provMat(sizMat, nbCol);
  std::vector<size_t> ListColSelect;
  std::vector<uint8_t> ListColSelect01(nbCol, 0);
  size_t eRank = 0;
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    f(provMat, eRank, iRow);
    for (size_t iRank = 0; iRank < eRank; iRank++) {
      size_t eCol = ListColSelect[iRank];
      T eVal1 = provMat(eRank, eCol);
      if (eVal1 != 0) {
        for (size_t iCol = eCol; iCol < nbCol; iCol++) {
          provMat(eRank, iCol) -= eVal1 * provMat(iRank, iCol);
        }
      }
    }
    auto get_firstnonzerocol_iife = [&]() -> size_t {
      for (size_t iCol = 0; iCol < nbCol; iCol++) {
        if (provMat(eRank, iCol) != 0)
          return iCol;
      }
      return std::numeric_limits<size_t>::max();
    };
    size_t FirstNonZeroCol = get_firstnonzerocol_iife();
    if (FirstNonZeroCol != std::numeric_limits<size_t>::max()) {
      ListColSelect.push_back(FirstNonZeroCol);
      ListColSelect01[FirstNonZeroCol] = 1;
      T eVal2 = 1 / provMat(eRank, FirstNonZeroCol);
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        provMat(eRank, iCol) *= eVal2;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        T eVal1 = provMat(iRank, FirstNonZeroCol);
        if (eVal1 != 0) {
          size_t StartCol = ListColSelect[iRank];
          for (size_t iCol = StartCol; iCol < nbCol; iCol++)
            provMat(iRank, iCol) -= eVal1 * provMat(eRank, iCol);
        }
      }
      eRank++;
    }
  }
  size_t nbVectZero = nbCol - eRank;
  MyMatrix<T> NSP = ZeroMatrix<T>(nbVectZero, nbCol);
  size_t nbVect = 0;
  for (size_t iCol = 0; iCol < nbCol; iCol++)
    if (ListColSelect01[iCol] == 0) {
      NSP(nbVect, iCol) = -1;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        size_t eCol = ListColSelect[iRank];
        NSP(nbVect, eCol) = provMat(iRank, iCol);
      }
      nbVect++;
    }
  return NSP;
}

//
template <typename T, typename F>
MyMatrix<T> NullspaceTrMatTarget_Kernel(size_t nbRow, size_t nbCol,
                                        size_t target_zero, F f) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in NullspaceTrMat_Kernel");
  size_t target_rank = nbCol - target_zero;
  MyMatrix<T> provMat(target_rank, nbCol);
  std::vector<size_t> ListColSelect;
  std::vector<uint8_t> ListColSelect01(nbCol, 0);
  size_t eRank = 0;
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    f(provMat, eRank, iRow);
    for (size_t iRank = 0; iRank < eRank; iRank++) {
      size_t eCol = ListColSelect[iRank];
      T eVal1 = provMat(eRank, eCol);
      if (eVal1 != 0) {
        for (size_t iCol = eCol; iCol < nbCol; iCol++) {
          provMat(eRank, iCol) -= eVal1 * provMat(iRank, iCol);
        }
      }
    }
    auto get_firstnonzerocol_iife = [&]() -> size_t {
      for (size_t iCol = 0; iCol < nbCol; iCol++) {
        if (provMat(eRank, iCol) != 0)
          return iCol;
      }
      return std::numeric_limits<size_t>::max();
    };
    size_t FirstNonZeroCol = get_firstnonzerocol_iife();
    if (FirstNonZeroCol != std::numeric_limits<size_t>::max()) {
      ListColSelect.push_back(FirstNonZeroCol);
      ListColSelect01[FirstNonZeroCol] = 1;
      T eVal2 = 1 / provMat(eRank, FirstNonZeroCol);
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        provMat(eRank, iCol) *= eVal2;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        T eVal1 = provMat(iRank, FirstNonZeroCol);
        if (eVal1 != 0) {
          size_t StartCol = ListColSelect[iRank];
          for (size_t iCol = StartCol; iCol < nbCol; iCol++)
            provMat(iRank, iCol) -= eVal1 * provMat(eRank, iCol);
        }
      }
      eRank++;
      if (eRank == target_rank) {
        // We reach the expected rank. So now returning.
        MyMatrix<T> NSP = ZeroMatrix<T>(target_zero, nbCol);
        size_t nbVect = 0;
        for (size_t iCol = 0; iCol < nbCol; iCol++) {
          if (ListColSelect01[iCol] == 0) {
            NSP(nbVect, iCol) = -1;
            for (size_t iRank = 0; iRank < eRank; iRank++) {
              size_t eCol = ListColSelect[iRank];
              NSP(nbVect, eCol) = provMat(iRank, iCol);
            }
            nbVect++;
          }
        }
        return NSP;
      }
    }
  }
  // The target was not achieved, we get a larger kernel than expected.
  // It could be a problem down the line, but such is life.
  size_t nbVectZero = nbCol - eRank;
  MyMatrix<T> NSP = ZeroMatrix<T>(nbVectZero, nbCol);
  size_t nbVect = 0;
  for (size_t iCol = 0; iCol < nbCol; iCol++)
    if (ListColSelect01[iCol] == 0) {
      NSP(nbVect, iCol) = -1;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        size_t eCol = ListColSelect[iRank];
        NSP(nbVect, eCol) = provMat(iRank, iCol);
      }
      nbVect++;
    }
  return NSP;
}

template <typename T, typename F>
MyVector<T> NullspaceTrMatTargetOne_Kernel(size_t nbRow, size_t nbCol, F f) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in NullspaceTrMat_Kernel");
  size_t target_rank = nbCol - 1;
  MyMatrix<T> provMat(target_rank, nbCol);
  std::vector<size_t> ListColSelect;
  std::vector<uint8_t> ListColSelect01(nbCol, 0);
  size_t eRank = 0;
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    f(provMat, eRank, iRow);
    for (size_t iRank = 0; iRank < eRank; iRank++) {
      size_t eCol = ListColSelect[iRank];
      T eVal1 = provMat(eRank, eCol);
      if (eVal1 != 0) {
        for (size_t iCol = eCol; iCol < nbCol; iCol++) {
          provMat(eRank, iCol) -= eVal1 * provMat(iRank, iCol);
        }
      }
    }
    auto get_firstnonzerocol_iife = [&]() -> size_t {
      for (size_t iCol = 0; iCol < nbCol; iCol++) {
        if (provMat(eRank, iCol) != 0)
          return iCol;
      }
      return std::numeric_limits<size_t>::max();
    };
    size_t FirstNonZeroCol = get_firstnonzerocol_iife();
    if (FirstNonZeroCol != std::numeric_limits<size_t>::max()) {
      ListColSelect.push_back(FirstNonZeroCol);
      ListColSelect01[FirstNonZeroCol] = 1;
      T eVal2 = 1 / provMat(eRank, FirstNonZeroCol);
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        provMat(eRank, iCol) *= eVal2;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        T eVal1 = provMat(iRank, FirstNonZeroCol);
        if (eVal1 != 0) {
          size_t StartCol = ListColSelect[iRank];
          for (size_t iCol = StartCol; iCol < nbCol; iCol++)
            provMat(iRank, iCol) -= eVal1 * provMat(eRank, iCol);
        }
      }
      eRank++;
      if (eRank == target_rank) {
        // We reach the expected rank. So now returning.
        MyVector<T> Vzero = ZeroVector<T>(nbCol);
        for (size_t iCol = 0; iCol < nbCol; iCol++) {
          if (ListColSelect01[iCol] == 0) {
            Vzero(iCol) = -1;
            for (size_t iRank = 0; iRank < eRank; iRank++) {
              size_t eCol = ListColSelect[iRank];
              Vzero(eCol) = provMat(iRank, iCol);
            }
            return Vzero;
          }
        }
      }
    }
  }
  // The target was not achieved, we get a larger kernel than expected.
  // We select one vector and it has to be processed down the line
  MyVector<T> Vzero = ZeroVector<T>(nbCol);
  auto set_vzero = [&]() -> void {
    for (size_t iCol = 0; iCol < nbCol; iCol++) {
      if (ListColSelect01[iCol] == 0) {
        Vzero(iCol) = -1;
        for (size_t iRank = 0; iRank < eRank; iRank++) {
          size_t eCol = ListColSelect[iRank];
          Vzero(eCol) = provMat(iRank, iCol);
        }
        return;
      }
    }
  };
  set_vzero();
  return Vzero;
}

/*
  This is for computing the kernel of a floating point matrix for which we know
  the rank. We want to avoid finding an incorrect kernel.
  ---
  Iniial equation is xM = 0
  Then we have xM M^T = 0
  Then we compute the lower eigenvalues and this gets us the approximate kernel.
 */
template <typename T>
MyMatrix<T> DominantKernel(MyMatrix<T> const &M, int const &k) {
  MyMatrix<T> Mprod = M * M.transpose();
  int dim_rel = Mprod.rows();
  Eigen::SelfAdjointEigenSolver<MyMatrix<T>> eig(Mprod);
  MyVector<T> ListEig = eig.eigenvalues();
  MyMatrix<T> ListVect = eig.eigenvectors();
#ifdef DEBUG_MAT_MATRIX
  for (int i = 1; i < dim_rel; i++) {
    if (ListEig(i - 1) > ListEig(i)) {
      std::cerr << "ListEig(i-1)=" << ListEig(i - 1)
                << " ListEig(i)=" << ListEig(i) << "\n";
      throw TerminalException{1};
    }
  }
#endif
  MyMatrix<T> NSP(k, dim_rel);
  for (int u = 0; u < k; u++) {
    for (int i = 0; i < dim_rel; i++) {
      NSP(u, i) = ListVect(i, u);
    }
  }
#ifdef DEBUG_MAT_MATRIX
  MyMatrix<T> prod = NSP * M;
  T error = L1_norm_mat(prod);
  if (error > 0.1) {
    std::cerr << "error=" << error << "\n";
    throw TerminalException{1};
  }
#endif
  return NSP;
}

template <typename T>
inline typename std::enable_if<is_ring_field<T>::value, MyMatrix<T>>::type
NullspaceTrMat(MyMatrix<T> const &Input) {
  size_t nbRow = Input.rows();
  size_t nbCol = Input.cols();
  auto f = [&](MyMatrix<T> &M, size_t eRank, size_t iRow) -> void {
    M.row(eRank) = Input.row(iRow);
  };
  return NullspaceTrMat_Kernel<T, decltype(f)>(nbRow, nbCol, f);
}

template <typename T>
inline typename std::enable_if<!is_ring_field<T>::value, MyMatrix<T>>::type
NullspaceTrMat(MyMatrix<T> const &Input) {
  // No division allowed. Maybe faster if we can allow for it using mpz_class.
  size_t nbRow = Input.rows();
  size_t nbCol = Input.cols();
  size_t maxRank = nbRow;
  if (nbCol < maxRank)
    maxRank = nbCol;
  size_t sizMat = maxRank + 1;
  MyMatrix<T> provMat(sizMat, nbCol);
  std::vector<size_t> ListColSelect;
  std::vector<uint8_t> ListColSelect01(nbCol, 0);
  size_t eRank = 0;
  size_t miss_val = std::numeric_limits<size_t>::max();
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    for (size_t iCol = 0; iCol < nbCol; iCol++)
      provMat(eRank, iCol) = Input(iRow, iCol);
    for (size_t iRank = 0; iRank < eRank; iRank++) {
      int eCol = ListColSelect[iRank];
      T eVal1 = provMat(eRank, eCol);
      T eVal2 = provMat(iRank, eCol);
      if (eVal1 != 0) {
        for (size_t iCol = 0; iCol < nbCol; iCol++)
          provMat(eRank, iCol) =
              provMat(eRank, iCol) * eVal2 - provMat(iRank, iCol) * eVal1;
      }
    }
    auto get_firstnonzerocol_iife = [&]() -> size_t {
      for (size_t iCol = 0; iCol < nbCol; iCol++)
        if (provMat(eRank, iCol) != 0)
          return iCol;
      return miss_val;
    };
    size_t FirstNonZeroCol = get_firstnonzerocol_iife();
    if (FirstNonZeroCol != miss_val) {
      ListColSelect.push_back(FirstNonZeroCol);
      ListColSelect01[size_t(FirstNonZeroCol)] = 1;
      T eVal2 = provMat(eRank, FirstNonZeroCol);
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        T eVal1 = provMat(iRank, FirstNonZeroCol);
        if (eVal1 != 0) {
          for (size_t iCol = 0; iCol < nbCol; iCol++)
            provMat(iRank, iCol) =
                provMat(iRank, iCol) * eVal2 - provMat(eRank, iCol) * eVal1;
        }
      }
      eRank++;
    }
  }
  // CODE INCOMPLETE BELOW. SOMETHING IS NEEDED.
  size_t nbVectZero = nbCol - eRank;
  MyMatrix<T> NSP = ZeroMatrix<T>(nbVectZero, nbCol);
  size_t nbVect = 0;
  for (size_t iCol = 0; iCol < nbCol; iCol++)
    if (ListColSelect01[iCol] == 0) {
      NSP(nbVect, iCol) = 1;
      T prodVal = 1;
      for (size_t iRank = 0; iRank < eRank; iRank++) {
        size_t eCol = ListColSelect[iRank];
        T pivotVal = provMat(iRank, iCol);
        if (pivotVal != 0) {
          for (size_t jRank = 0; jRank < iRank; jRank++) {
            size_t fCol = ListColSelect[jRank];
            NSP(nbVect, fCol) *= provMat(iRank, eCol);
          }
          NSP(nbVect, iCol) *= provMat(iRank, eCol);
          //
          NSP(nbVect, eCol) = -prodVal * provMat(iRank, iCol);
          prodVal *= provMat(iRank, eCol);
        }
      }
      nbVect++;
    }
  return NSP;
}

template <typename T> MyMatrix<T> NullspaceMat(MyMatrix<T> const &M) {
  return TMat_SelectRowCol(TransposedMat(M)).NSP;
}

template <typename T> int RankMatKernel(MyMatrix<T> const &Input) {
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(Input);
  return eSelect.TheRank;
}

template <typename T>
inline typename std::enable_if<is_ring_field<T>::value, int>::type
RankMat(MyMatrix<T> const &Input) {
  return RankMatKernel(Input);
}

template <typename T>
inline typename std::enable_if<!is_ring_field<T>::value, int>::type
RankMat(MyMatrix<T> const &Input) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> InputF = UniversalMatrixConversion<Tfield, T>(Input);
  return RankMatKernel(InputF);
}

template <typename T> MyMatrix<T> IdentityMat(int n) {
  T t;
  MyMatrix<T> TheMat(n, n);
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      if (i == j)
        t = 1;
      else
        t = 0;
      TheMat(i, j) = t;
    }
  return TheMat;
}

template <typename T>
void TMat_ImageIntVector(MyVector<T> &eVect, MyMatrix<T> &TheMat,
                         MyVector<T> &eVectImg) {
  int iCol, nbCol, iRow, nbRow, n;
  n = eVect->n;
  nbRow = TheMat.rows();
  nbCol = TheMat.cols();
  if (n != nbRow) {
    std::cerr << "Error in ImageIntVector\n";
    std::cerr << "n=" << n << " nbRow=" << nbRow << "\n";
    throw TerminalException{1};
  }
  for (iCol = 0; iCol < nbCol; iCol++) {
    T t = 0;
    for (iRow = 0; iRow < n; iRow++)
      t += TheMat(iRow, iCol) * eVect(iRow);
    eVectImg(iCol) = t;
  }
}

template <typename T> MyMatrix<T> CongrMap(MyMatrix<T> const &eMat) {
  MyMatrix<T> TheInv = Inverse(eMat);
  return TransposedMat(TheInv);
}

template <typename T> T DeterminantMatKernel(MyMatrix<T> const &TheMat) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in DeterminantMatKernel");
  T hVal, alpha;
  int n = TheMat.rows();
  MyMatrix<T> WorkMat = TheMat;
  std::vector<int> eVectPos(n, -1);
  T TheDet(1);
  for (int i = 0; i < n; i++) {
    int jSel = -1;
    for (int j = 0; j < n; j++)
      if (eVectPos[j] == -1) {
        hVal = WorkMat(i, j);
        if (hVal != 0)
          jSel = j;
      }
    if (jSel == -1) {
      return T(0);
    }
    eVectPos[jSel] = i;
    for (int j = 0; j < n; j++)
      if (j != jSel) {
        alpha = WorkMat(i, j) / WorkMat(i, jSel);
        for (int k = 0; k < n; k++)
          WorkMat(k, j) -= alpha * WorkMat(k, jSel);
      }
    TheDet = TheDet * WorkMat(i, jSel);
  }
  int nbchg = 0;
  for (int i = 0; i < n - 1; i++)
    for (int j = i + 1; j < n; j++)
      if (eVectPos[i] > eVectPos[j])
        nbchg++;
  int res = nbchg % 2;
  if (res == 0)
    return TheDet;
  return -TheDet;
}

template <typename T>
inline typename std::enable_if<is_ring_field<T>::value, T>::type
DeterminantMat(MyMatrix<T> const &Input) {
  return DeterminantMatKernel(Input);
}

template <typename T>
inline typename std::enable_if<!is_ring_field<T>::value, T>::type
DeterminantMat(MyMatrix<T> const &Input) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> InputF = UniversalMatrixConversion<Tfield, T>(Input);
  Tfield eDet_field = DeterminantMatKernel(InputF);
  return UniversalScalarConversion<T, Tfield>(eDet_field);
}

// A significantly slower algorithm for computing the determinant.
// It is good as a control for the above method and can also be used
// for consistency checks of arithmetics.
template <typename T> T DeterminantMatPermutation(MyMatrix<T> const &A) {
  int n = A.rows();
  if (n == 0)
    return T(1);
  std::vector<int> s(n);
  for (int i = 0; i < n; i++)
    s[i] = i;
  T TheDet(0);
  do {
    T eProd(1);
    for (int u = 0; u < n; u++)
      eProd *= A(u, s[u]);
    int eSign = 1;
    for (int i = 0; i < n; i++)
      for (int j = i + 1; j < n; j++)
        if (s[j] < s[i])
          eSign = -eSign;
    TheDet += eSign * eProd;
  } while (std::next_permutation(s.begin(), s.end()));
  return TheDet;
}

template <typename T> MyVector<T> OrthogonalHyperplane(MyMatrix<T> const &M) {
  int nbVert = M.rows();
  MyVector<T> eVect(nbVert);
  MyMatrix<T> Mred(nbVert - 1, nbVert - 1);
  int eCoeff = 1;
  for (int iVert = 0; iVert < nbVert; iVert++) {
    int iRow = 0;
    for (int iLine = 0; iLine < nbVert; iLine++)
      if (iLine != iVert) {
        for (int iCol = 0; iCol < nbVert - 1; iCol++)
          Mred(iRow, iCol) = M(iLine, iCol);
        iRow++;
      }
    eVect(iVert) = eCoeff * DeterminantMat(Mred);
    eCoeff *= -1;
  }
  return eVect;
}

template <typename T>
MyMatrix<T> SelectRow(MyMatrix<T> const &TheMat,
                      std::vector<int> const &eList) {
  size_t nbRowRed = eList.size();
  size_t nbCol = TheMat.cols();
  MyMatrix<T> TheProv(nbRowRed, nbCol);
  for (size_t iRow = 0; iRow < nbRowRed; iRow++) {
    size_t jRow = eList[iRow];
    TheProv.row(iRow) = TheMat.row(jRow);
  }
  return TheProv;
}

template <typename T>
MyMatrix<T> SelectColumn(MyMatrix<T> const &TheMat,
                         std::vector<int> const &eList) {
  size_t nbRow = TheMat.rows();
  size_t nbColRed = eList.size();
  MyMatrix<T> TheProv(nbRow, nbColRed);
  for (size_t iCol = 0; iCol < nbColRed; iCol++) {
    size_t jCol = eList[iCol];
    TheProv.col(iCol) = TheMat.col(jCol);
  }
  return TheProv;
}

template <typename T>
MyVector<T> SelectColumnVector(MyVector<T> const &TheV,
                               std::vector<int> const &eList) {
  int nbColRed = eList.size();
  MyVector<T> TheProv(nbColRed);
  for (int iCol = 0; iCol < nbColRed; iCol++) {
    int jCol = eList[iCol];
    TheProv(iCol) = TheV(jCol);
  }
  return TheProv;
}

template <typename T>
bool TestEquality(MyVector<T> const &V1, MyVector<T> const &V2) {
  int n1 = V1.size();
  int n2 = V2.size();
  if (n1 != n2)
    return false;
  for (int i = 0; i < n1; i++)
    if (V1(i) != V2(i))
      return false;
  return true;
}

template <typename T>
bool TestEqualityMatrix(MyMatrix<T> const &M1, MyMatrix<T> const &M2) {
  int n1 = M1.size();
  int n2 = M2.size();
  if (n1 != n2)
    return false;
  for (int i = 0; i < n1; i++)
    if (M1(i) != M2(i))
      return false;
  return true;
}

// In input a matrix
// in output the lines that spann its rank in sequential order
template <typename T> MyMatrix<T> RowReduction(MyMatrix<T> const &eMatIn) {
  SelectionRowCol<T> eSelect = TMat_SelectRowCol(eMatIn);
  std::vector<int> ListRowSelect = eSelect.ListRowSelect;
  return SelectRow(eMatIn, ListRowSelect);
}

template <typename T> bool IsZeroVector(MyVector<T> const &V) {
  int n = V.size();
  for (int i = 0; i < n; i++)
    if (V(i) != 0)
      return false;
  return true;
}

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

template <typename T> bool IsIntegerVector(MyVector<T> const &V) {
  int n = V.size();
  for (int i = 0; i < n; i++)
    if (!IsInteger(V(i)))
      return false;
  return true;
}

template <typename T>
inline typename std::enable_if<is_ring_field<T>::value,
                               std::optional<MyVector<T>>>::type
SolutionMat(MyMatrix<T> const &eMat, MyVector<T> const &eVect) {
  return SolutionMatKernel(eMat, eVect);
}

template <typename T>
inline typename std::enable_if<!is_ring_field<T>::value,
                               std::optional<MyVector<T>>>::type
SolutionMat(MyMatrix<T> const &eMat, MyVector<T> const &eVect) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> eMatF = UniversalMatrixConversion<Tfield, T>(eMat);
  MyVector<Tfield> eVectF = UniversalVectorConversion<Tfield, T>(eVect);
  std::optional<MyVector<Tfield>> opt = SolutionMatKernel(eMatF, eVectF);
  if (opt) {
    const MyVector<Tfield> &V = *opt;
    if (!IsIntegerVector(V))
      return {};
    return UniversalVectorConversion<T, Tfield>(V);
  }
  return {};
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

template <typename T> struct MutableSubspaceBelongingRepetitive {
private:
  size_t dim;
  size_t n_vect;
  std::vector<T> TheBasis;
  std::vector<size_t> ListColSelect;
  std::vector<T> Vins;

public:
  MutableSubspaceBelongingRepetitive(int const &_dim)
      : dim(_dim), n_vect(0), Vins(dim) {}
  void SetVectout(int i_vect, MyVector<T> &Vout) {
    size_t pos = i_vect * dim;
    for (size_t i = 0; i < dim; i++) {
      Vout(i) = TheBasis[pos];
      pos++;
    }
  }
  MyMatrix<T> GetBasis() const {
    MyMatrix<T> Mret(n_vect, dim);
    size_t pos = 0;
    for (size_t i_vect = 0; i_vect < n_vect; i_vect++) {
      for (size_t u = 0; u < dim; u++) {
        Mret(i_vect, u) = TheBasis[pos];
        pos++;
      }
    }
    return Mret;
  }
  size_t get_n_vect() const { return n_vect; }

private:
  size_t IsInspace_priv() {
    size_t pos = 0;
    for (size_t i_vect = 0; i_vect < n_vect; i_vect++) {
      size_t eCol = ListColSelect[i_vect];
      T coeff = Vins[eCol];
      if (coeff != 0) {
        for (size_t i = 0; i < dim; i++) {
          Vins[i] -= coeff * TheBasis[pos];
          pos++;
        }
      } else {
        pos += dim;
      }
    }
    for (size_t i = 0; i < dim; i++) {
      if (Vins[i] != 0) {
        return i;
      }
    }
    return dim;
  }

public:
  bool IsInspace(MyVector<T> const &v) {
    for (size_t i = 0; i < dim; i++) {
      Vins[i] = v(i);
    }
    return IsInspace_priv() == dim;
  }

private:
  bool InsertVector_priv() {
    size_t idx = IsInspace_priv();
    if (idx == dim) {
      return true;
    }
    T coeff = Vins[idx];
    for (size_t i = 0; i < dim; i++) {
      Vins[i] /= coeff;
    }
    for (size_t i_vect = 0; i_vect < n_vect; i_vect++) {
      T coeff = TheBasis[i_vect * dim + idx];
      if (coeff != 0) {
        for (size_t i = 0; i < dim; i++) {
          TheBasis[i_vect * dim + i] -= coeff * Vins[i];
        }
      }
    }
    ListColSelect.push_back(idx);
    n_vect++;
    for (size_t u = 0; u < dim; u++) {
      TheBasis.push_back(Vins[u]);
    }
    return false;
  }

public:
  bool InsertVector(MyVector<T> const &v) {
    for (size_t i = 0; i < dim; i++) {
      Vins[i] = v(i);
    }
    return InsertVector_priv();
  }
  void InsertMatrix(MyMatrix<T> const &M) {
    int n_row = M.rows();
    for (int i_row = 0; i_row < n_row; i_row++) {
      for (size_t i = 0; i < dim; i++) {
        Vins[i] = M(i_row, i);
      }
      (void)InsertVector_priv();
    }
  }
};

// Returns true if M1 is a subspace of M2
template <typename T>
bool IsSubspaceContained(MyMatrix<T> const &M1, MyMatrix<T> const &M2) {
  int dim = M2.cols();
  MutableSubspaceBelongingRepetitive<T> msbr(dim);
  msbr.InsertMatrix(M2);
  for (int iRow = 0; iRow < M1.rows(); iRow++) {
    MyVector<T> V = GetMatrixRow(M1, iRow);
    if (!msbr.IsInspace(V)) {
      return false;
    }
  }
  return true;
}

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
MyVector<T> GetMatrixRow(MyMatrix<T> const &M, int const &iRow) {
  int nbCol = M.cols();
  MyVector<T> V(nbCol);
  for (int iCol = 0; iCol < nbCol; iCol++)
    V(iCol) = M(iRow, iCol);
  return V;
}

template <typename T> MyMatrix<T> MatrixFromVector(MyVector<T> const &V) {
  int n = V.size();
  MyMatrix<T> M(1, n);
  for (int i = 0; i < n; i++) {
    M(0, i) = V(i);
  }
  return M;
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

template <typename T> MyMatrix<T> SelectNonZeroRows(MyMatrix<T> const &EXT) {
  int nbRow = EXT.rows();
  int nbCol = EXT.cols();
  auto IsLineZero = [&](int const &iLine) -> bool {
    for (int iCol = 0; iCol < nbCol; iCol++)
      if (EXT(iLine, iCol) != 0)
        return false;
    return true;
  };
  std::vector<int> ListIdx;
  for (int iRow = 0; iRow < nbRow; iRow++) {
    if (!IsLineZero(iRow))
      ListIdx.push_back(iRow);
  }
  return SelectRow(EXT, ListIdx);
}

template <typename T>
std::vector<int> ColumnReductionSet_Kernel(MyMatrix<T> const &eMatIn) {
  std::vector<int> ListCol = TMat_SelectRowCol(eMatIn).ListColSelect;
  std::sort(ListCol.begin(), ListCol.end());
  return ListCol;
}

template <typename T>
inline typename std::enable_if<is_ring_field<T>::value, std::vector<int>>::type
ColumnReductionSet(MyMatrix<T> const &Input) {
  return ColumnReductionSet_Kernel(Input);
}

template <typename T>
inline typename std::enable_if<!is_ring_field<T>::value, std::vector<int>>::type
ColumnReductionSet(MyMatrix<T> const &Input) {
  using Tfield = typename overlying_field<T>::field_type;
  MyMatrix<Tfield> InputF = UniversalMatrixConversion<Tfield, T>(Input);
  return ColumnReductionSet_Kernel(InputF);
}

template <typename T>
MyMatrix<T> ColRowSymmetricMatrix(MyMatrix<T> const &M,
                                  std::vector<int> const &LSel) {
  size_t siz = LSel.size();
  MyMatrix<T> Mred(siz, siz);
  for (size_t j = 0; j < siz; j++)
    for (size_t i = 0; i < siz; i++)
      Mred(i, j) = M(LSel[i], LSel[j]);
  return Mred;
}

template <typename T> MyMatrix<T> ColumnReduction(MyMatrix<T> const &eMatIn) {
  std::vector<int> l_cols = ColumnReductionSet(eMatIn);
  return SelectColumn(eMatIn, l_cols);
}

template <typename T> MyMatrix<T> ExtendToBasis(MyMatrix<T> const &M) {
  int n_cols = M.cols();
  int n_rows = M.rows();
  std::vector<int> V = TMat_SelectRowCol(M).ListColSelect;
  if (size_t(n_rows) != V.size()) {
    std::cerr << "M=\n";
    WriteMatrix(std::cerr, M);
    std::cerr << "rank=" << RankMat(M) << " n_cols=" << n_cols
              << " n_row=" << n_rows << "\n";
    std::cerr << "V =";
    for (auto &val : V)
      std::cerr << " " << val;
    std::cerr << "\n";
    std::cerr
        << "The original matrix M does not appear to be linearly independent\n";
    throw TerminalException{1};
  }
  MyMatrix<T> Mret(n_cols, n_cols);
  for (int i_row = 0; i_row < n_rows; i_row++)
    Mret.row(i_row) = M.row(i_row);
  std::vector<uint8_t> f(n_cols, 0);
  for (auto &val : V)
    f[val] = 1;
  size_t pos = n_rows;
  for (int i_col = 0; i_col < n_cols; i_col++) {
    if (f[i_col] == 0) {
      for (int j_col = 0; j_col < n_cols; j_col++) {
        if (i_col == j_col)
          Mret(pos, j_col) = 1;
        else
          Mret(pos, j_col) = 0;
      }
      pos++;
    }
  }
  return Mret;
}

template <typename T>
MyMatrix<T> Concatenate(MyMatrix<T> const &eMat1, MyMatrix<T> const &eMat2) {
  int nbRow1 = eMat1.rows();
  int nbRow2 = eMat2.rows();
  int nbCol1 = eMat1.cols();
  int nbCol2 = eMat2.cols();
  if (nbCol1 != nbCol2) {
    std::cerr << "nbCol1=" << nbCol1 << " nbCol2=" << nbCol2 << "\n";
    std::cerr << "Error in the number of columns\n";
    throw TerminalException{1};
  }
  int nbCol = nbCol1;
  MyMatrix<T> eMatRet(nbRow1 + nbRow2, nbCol);
  for (int iCol = 0; iCol < nbCol; iCol++)
    for (int iRow = 0; iRow < nbRow1; iRow++)
      eMatRet(iRow, iCol) = eMat1(iRow, iCol);
  for (int iCol = 0; iCol < nbCol; iCol++)
    for (int iRow = 0; iRow < nbRow2; iRow++)
      eMatRet(iRow + nbRow1, iCol) = eMat2(iRow, iCol);
  return eMatRet;
}

template <typename T>
MyMatrix<T> ConcatenateMatVec(MyMatrix<T> const &M, MyVector<T> const &V) {
  int nbRow = M.rows();
  int nbColM = M.cols();
  int nbCol = V.size();
  if (nbRow != 0) {
    if (nbColM != nbCol) {
      std::cerr << "Error in ConcatenateMatVec\n";
      std::cerr << "We have nbCol=" << nbCol << " nbColM=" << nbColM << "\n";
      throw TerminalException{1};
    }
  }
  MyMatrix<T> Mret(nbRow + 1, nbCol);
  for (int iCol = 0; iCol < nbCol; iCol++)
    for (int iRow = 0; iRow < nbRow; iRow++)
      Mret(iRow, iCol) = M(iRow, iCol);
  for (int iCol = 0; iCol < nbCol; iCol++)
    Mret(nbRow, iCol) = V(iCol);
  return Mret;
}

template <typename T>
MyMatrix<T> ConcatenateMatVec_Tr(MyMatrix<T> const &M, MyVector<T> const &V) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  int n = V.size();
  if (nbCol != 0) {
    if (nbRow != n) {
      std::cerr << "Error in ConcatenateMatVec_Tr\n";
      std::cerr << "We have nbRow=" << nbRow << " n=" << n << "\n";
      throw TerminalException{1};
    }
  }
  MyMatrix<T> Mret(nbRow, nbCol + 1);
  for (int iCol = 0; iCol < nbCol; iCol++)
    for (int iRow = 0; iRow < nbRow; iRow++)
      Mret(iRow, iCol) = M(iRow, iCol);
  for (int i = 0; i < n; i++)
    Mret(i, nbCol) = V(i);
  return Mret;
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
  MyMatrix<T> M = ZeroMatrix<T>(nbRow, nbCol);
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

template <typename T>
T MatrixScalarProduct(MyMatrix<T> const &M1, MyMatrix<T> const &M2) {
  int n1 = M1.rows();
  int p1 = M1.cols();
  int n2 = M2.rows();
  int p2 = M2.cols();
  if (n1 != n2 || p1 != p2) {
    std::cerr << "Incoherency\n";
    std::cerr << "n1=" << n1 << " n2=" << n2 << "\n";
    std::cerr << "p1=" << p1 << " p1=" << p2 << "\n";
    throw TerminalException{1};
  }
  T eSum = 0;
  for (int i = 0; i < n1; i++)
    for (int j = 0; j < p1; j++)
      eSum += M1(i, j) * M2(i, j);
  return eSum;
}

template <typename T>
MyMatrix<T> VectorToSymmetricMatrix(MyVector<T> const &V, int const &n) {
  MyMatrix<T> RetMat = ZeroMatrix<T>(n, n);
  int idx = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++) {
      RetMat(i, j) = V(idx);
      RetMat(j, i) = V(idx);
      idx++;
    }
  return RetMat;
}

template <typename T>
MyMatrix<T> VectorToSymmetricMatrixB(MyVector<T> const &V, int const &n) {
  MyMatrix<T> RetMat(n, n);
  int idx = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++) {
      T eVal = V(idx);
      if (i != j)
        eVal = eVal / 2;
      RetMat(i, j) = eVal;
      RetMat(j, i) = eVal;
      idx++;
    }
  return RetMat;
}

template <typename T>
MyVector<T> SymmetricMatrixToVector(MyMatrix<T> const &M) {
  int n = M.rows();
  int dim = (n * (n + 1)) / 2;
  MyVector<T> eVect(dim);
  int idx = 0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      eVect(idx) = M(i, j);
      idx++;
    }
  }
  return eVect;
}

template<typename T>
MyVector<T> MatrixToVector(MyMatrix<T> const& M) {
  int n = M.rows();
  int dim = n * n;
  MyVector<T> V(dim);
  int idx = 0;
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      V(idx) = M(i,j);
      idx += 1;
    }
  }
  return V;
}

template<typename T>
MyMatrix<T> VectorToMatrix(MyVector<T> const& V, int const& n) {
  MyMatrix<T> M(n, n);
  int idx=0;
  for (int i=0; i<n; i++) {
    for (int j=0; j<n; j++) {
      M(i,j) = V(idx);
      idx += 1;
    }
  }
  return M;
}


template <typename T> MyVector<T> GetSymmetricMatrixWeightVector(int const &n) {
  int dim = (n * (n + 1)) / 2;
  MyVector<T> eVect(dim);
  int idx = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++) {
      T eVal;
      if (i == j)
        eVal = 1;
      else
        eVal = 2;
      eVect(idx) = eVal;
      idx++;
    }
  return eVect;
}

template <typename T>
MyVector<T> SymmetricMatrixToVectorB(MyMatrix<T> const &M) {
  int n = M.rows();
  int dim = (n * (n + 1)) / 2;
  MyVector<T> eVect(dim);
  int idx = 0;
  T eVal;
  for (int i = 0; i < n; i++)
    for (int j = 0; j <= i; j++) {
      eVal = M(i, j);
      if (i != j)
        eVal = 2 * eVal;
      eVect(idx) = eVal;
    }
  return eVect;
}

template <typename T>
void PrintEigenvalueDefect(MyMatrix<T> const &Sinp,
                           MyVector<T> const &ListEigVal,
                           MyMatrix<T> const &ListEigVect, std::ostream &os) {
  int n = Sinp.rows();
  for (int i = 0; i < n; i++) {
    T eDelta = 0;
    for (int j = 0; j < n; j++) {
      T eSum = ListEigVal(i) * ListEigVect(i, j);
      for (int k = 0; k < n; k++)
        eSum -= ListEigVect(i, k) * Sinp(j, k);
      eDelta += T_abs(eSum);
    }
    os << "i=" << i << " err=" << eDelta << "\n";
  }
}

template <typename T>
MyVector<T> SolveConjGrad(MyMatrix<T> const &A, MyVector<T> const &b) {
  static_assert(is_ring_field<T>::value,
                "Requires T to be a field in SolveConjGrad");
  int n = b.size();
  MyVector<T> r(n);
  MyVector<T> p(n);
  MyVector<T> x(n);
  MyVector<T> Ap(n);
  double rsold, rsnew, alpha;
  r = b;
  p = r;
  rsold = r.dot(r);
  T eZer = 0;
  for (int i = 0; i < n; i++)
    x(i) = eZer;
  int nbOper = 4 * n;
  for (int i = 0; i < nbOper; i++) {
    Ap = A * p;
    alpha = rsold / EvaluationQuadForm(A, p);
    x += alpha * p;
    r -= alpha * Ap;
    rsnew = r.dot(r);
    p = r + (rsnew / rsold) * p;
    rsold = rsnew;
  }
  return x;
}

template <typename T>
MyVector<T> Concatenation(MyVector<T> const &V1, MyVector<T> const &V2) {
  int dim1 = V1.size();
  int dim2 = V2.size();
  MyVector<T> V(dim1 + dim2);
  for (int i = 0; i < dim1; i++)
    V(i) = V1(i);
  for (int i = 0; i < dim2; i++)
    V(i + dim1) = V2(i);
  return V;
}

template <typename T> MyVector<T> Isobarycenter(MyMatrix<T> const &eMat) {
  int nbRow = eMat.rows();
  int nbCol = eMat.cols();
  MyVector<T> eVect(nbCol);
  T nbRow_T = nbRow;
  for (int iCol = 0; iCol < nbCol; iCol++) {
    T eSum = 0;
    for (int iRow = 0; iRow < nbRow; iRow++)
      eSum += eMat(iRow, iCol);
    eVect(iCol) = eSum / nbRow_T;
  }
  return eVect;
}

template <typename T> bool IsZeroMatrix(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  for (int iCol = 0; iCol < nbCol; iCol++)
    for (int iRow = 0; iRow < nbRow; iRow++)
      if (M(iRow, iCol) != 0)
        return false;
  return true;
}

template <typename T>
inline typename std::enable_if<is_ring_field<T>::value, MyVector<T>>::type
CanonicalizeVector(MyVector<T> const &V) {
  int n = V.size();
  T TheMin(0);
  int iSelect = -1;
  for (int i = 0; i < n; i++) {
    T eVal = V(i);
    if (eVal != 0) {
      T eAbs = T_abs(eVal);
      if (iSelect == -1) {
        TheMin = eAbs;
        iSelect = i;
      } else {
        if (eAbs < TheMin) {
          TheMin = eAbs;
          iSelect = i;
        }
      }
    }
  }
  if (iSelect == -1)
    return V;
  T eQuot = 1 / TheMin;
  return eQuot * V;
}

/* return true if V1 < V2 according to lexicographic order */
template <typename T>
bool IsLower(MyVector<T> const &V1, MyVector<T> const &V2) {
  int n = V1.size();
  for (int i = 0; i < n; i++) {
    if (V1(i) != V2(i)) {
      if (V1(i) < V2(i)) {
        return true;
      } else {
        return false;
      }
    }
  }
  return false;
}

template <typename T>
bool operator==(MyVector<T> const &V1, MyVector<T> const &V2) {
  int n = V1.size();
  if (V2.size() != n) {
    std::cerr << "We should not have different sizes\n";
    throw TerminalException{1};
  }
  for (int i = 0; i < n; i++) {
    if (V1(i) != V2(i))
      return false;
  }
  return true;
}

template <typename T>
bool operator<(MyMatrix<T> const &M1, MyMatrix<T> const &M2) {
  int nbRow = M1.rows();
  int nbCol = M1.cols();
  for (int iRow = 0; iRow < nbRow; iRow++)
    for (int iCol = 0; iCol < nbCol; iCol++) {
      if (M1(iRow, iCol) < M2(iRow, iCol))
        return true;
      if (M1(iRow, iCol) > M2(iRow, iCol))
        return false;
    }
  return false;
}

template <typename T>
bool operator<(MyVector<T> const &V1, MyVector<T> const &V2) {
  int siz = V1.size();
  for (int i = 0; i < siz; i++) {
    if (V1(i) < V2(i))
      return true;
    if (V1(i) > V2(i))
      return false;
  }
  return false;
}

namespace std {
template <typename T> struct less<MyVector<T>> {
  bool operator()(MyVector<T> const &V1, MyVector<T> const &V2) const {
    int siz = V1.size();
    for (int i = 0; i < siz; i++) {
      if (V1(i) < V2(i))
        return true;
      if (V2(i) < V1(i))
        return false;
    }
    return false;
  }
};
// clang-format off
}  // namespace std
// clang-format on

template <typename T> T L1_norm_vect(MyVector<T> const &V) {
  int siz = V.size();
  T norm = 0;
  for (int i = 0; i < siz; i++)
    norm += T_abs(V(i));
  return norm;
}

template <typename T> T Linfinity_norm_vect(MyVector<T> const &V) {
  int siz = V.size();
  T norm = 0;
  for (int i = 0; i < siz; i++)
    norm = T_max(norm, T_abs(V(i)));
  return norm;
}

template <typename T> T L1_norm_mat(MyMatrix<T> const &M) {
  int siz = M.size();
  T norm = 0;
  for (int i = 0; i < siz; i++)
    norm += T_abs(M(i));
  return norm;
}

template <typename T> T Linfinity_norm_mat(MyMatrix<T> const &M) {
  int siz = M.size();
  T norm = 0;
  for (int i = 0; i < siz; i++)
    norm = T_max(norm, T_abs(M(i)));
  return norm;
}

// Just sorting (no unicization)
//
template <typename T> MyMatrix<T> SortMatrix(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  auto comp = [](MyVector<T> const &V1, MyVector<T> const &V2) -> bool {
    return IsLower(V1, V2);
  };
  std::set<MyVector<T>, std::function<bool(MyVector<T>, MyVector<T>)>>
      eListVect(comp);
  //  eListVect=std::set<MyVector<T>, decltype(comp)> (comp);
  for (int iRow = 0; iRow < nbRow; iRow++) {
    MyVector<T> V(nbCol);
    for (int iCol = 0; iCol < nbCol; iCol++)
      V(iCol) = M(iRow, iCol);
    eListVect.insert(V);
  }
  int nbVect = eListVect.size();
  MyMatrix<T> Mret(nbVect, nbCol);
  int idx = 0;
  for (auto &eV : eListVect) {
    for (int iCol = 0; iCol < nbCol; iCol++)
      Mret(idx, iCol) = eV(iCol);
    idx++;
  }
  return Mret;
}

template <typename T>
MyVector<T> GetMatrixCol(MyMatrix<T> const &M, int const &iCol) {
  int nbRow = M.rows();
  MyVector<T> V(nbRow);
  for (int iRow = 0; iRow < nbRow; iRow++)
    V(iRow) = M(iRow, iCol);
  return V;
}

template <typename T>
std::vector<T> ConcatenateVect(std::vector<T> const &ListV1,
                               std::vector<T> const &ListV2) {
  std::vector<T> ListV = ListV1;
  for (auto &eV : ListV2)
    ListV.push_back(eV);
  return ListV;
}

template <typename T>
MyMatrix<T> MatrixFromVectorFamily(std::vector<MyVector<T>> const &ListVect) {
  int nbVect = ListVect.size();
  if (nbVect == 0) {
    std::cerr << "Error in MatrixFromVectorFamily. ListVect is empty\n";
    std::cerr
        << "We cannot create the matrix since we cannot know the dimension\n";
    throw TerminalException{1};
  }
  size_t dim = ListVect[0].size();
  int dim_i = static_cast<int>(dim);
  MyMatrix<T> M(nbVect, dim_i);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    if (ListVect[iVect].size() != dim_i) {
      std::cerr << "Vector lengths are not homogeneous\n";
      throw TerminalException{1};
    }
    for (int i = 0; i < dim_i; i++)
      M(iVect, i) = ListVect[iVect](i);
  }
  return M;
}

template<typename T>
MyMatrix<T> MatrixFromPairVector(std::pair<MyVector<T>, MyVector<T>> const& pair) {
  int n1 = pair.first.size();
  int n2 = pair.second.size();
  if (n1 != n2) {
    std::cerr << "n1=" << n1 << " n2=" << n2 << "\n";
    std::cerr << "Vector lengths are not homogeneous\n";
    throw TerminalException{1};
  }
  int n = n1;
  MyMatrix<T> M(2,n);
  for (int i=0; i<n; i++) {
    M(0, i) = pair.first(i);
    M(1, i) = pair.second(i);
  }
  return M;
}


template <typename T>
MyMatrix<T>
MatrixFromVectorFamilyDim(int const &dim,
                          std::vector<MyVector<T>> const &ListVect) {
  int nbVect = ListVect.size();
  MyMatrix<T> M(nbVect, dim);
  for (int iVect = 0; iVect < nbVect; iVect++) {
    for (int i = 0; i < dim; i++)
      M(iVect, i) = ListVect[iVect](i);
  }
  return M;
}

template <typename T> MyVector<T> SumMatrix(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  MyVector<T> V(nbCol);
  for (int iCol = 0; iCol < nbCol; iCol++) {
    T eSum = 0;
    for (int iRow = 0; iRow < nbRow; iRow++)
      eSum += M(iRow, iCol);
    V(iCol) = eSum;
  }
  return V;
}

template <typename T>
std::vector<T> StdVectorFromVector(MyVector<T> const &eV) {
  int siz = eV.size();
  std::vector<T> eVect(siz);
  for (int i = 0; i < siz; i++)
    eVect[i] = eV(i);
  return eVect;
}

template <typename T>
MyVector<T> VectorFromStdVector(std::vector<T> const &eList) {
  int siz = eList.size();
  MyVector<T> eVect(siz);
  for (int i = 0; i < siz; i++)
    eVect(i) = eList[i];
  return eVect;
}

template <typename T>
inline typename std::enable_if<is_float_arithmetic<T>::value, MyVector<T>>::type
RemoveFractionVector(MyVector<T> const &V) {
  return V;
}

template <typename T, typename Tint>
T EvaluationQuadForm(MyMatrix<T> const &eMat, MyVector<Tint> const &eVect) {
  size_t n = eVect.size();
#ifdef SANITY_CHECK_MAT_MATRIX
  if (eMat.rows() != eVect.size() || eMat.cols() != eVect.size()) {
    std::cerr << "MAT: We should have eMat a square matrix of size n x n with n the length of eVect\n";
    throw TerminalException{1};
  }
#endif
  T eSum(0);
  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
      eSum += eVect(i) * eVect(j) * eMat(i, j);
  return eSum;
}

//
// The matrix is supposed to describe a vector space. We want
// to canonicalize it in order to get a better representation
// of the space.
//
template <typename T>
MyMatrix<T> CanonicalizeBasisVectorSpace(MyMatrix<T> const &inputMat) {
  MyMatrix<T> WorkMat = inputMat;
  size_t nbRow = WorkMat.rows();
  size_t nbCol = WorkMat.cols();
  std::vector<int> ColStatus(nbCol, 1);
  std::vector<int> RowStatus(nbRow, 1);
  for (size_t iRow = 0; iRow < nbRow; iRow++) {
    size_t FoundRow = std::numeric_limits<size_t>::max();
    size_t FoundCol = 0;
    T MaxValue = 0;
    for (size_t eRow = 0; eRow < nbRow; eRow++) {
      for (size_t eCol = 0; eCol < nbCol; eCol++) {
        if (ColStatus[eCol] == 1 && RowStatus[eRow] == 1) {
          T aVal = T_abs(WorkMat(eRow, eCol));
          if (FoundRow == std::numeric_limits<size_t>::max()) {
            MaxValue = aVal;
            FoundRow = eRow;
            FoundCol = eCol;
          } else {
            if (aVal > MaxValue) {
              MaxValue = aVal;
              FoundRow = eRow;
              FoundCol = eCol;
            }
          }
        }
      }
    }
    //
    ColStatus[FoundCol] = 0;
    RowStatus[FoundRow] = 0;
    for (size_t eRow = 0; eRow < nbRow; eRow++) {
      if (eRow != FoundRow) {
        T alpha = WorkMat(eRow, FoundCol) / WorkMat(FoundRow, FoundCol);
        for (size_t iCol = 0; iCol < nbCol; iCol++)
          WorkMat(eRow, iCol) -= alpha * WorkMat(FoundRow, iCol);
      }
    }
  }
  return WorkMat;
}

template <typename T, typename Tint>
T ScalarProductQuadForm(MyMatrix<T> const &eMat, MyVector<Tint> const &V1,
                        MyVector<Tint> const &V2) {
  int n = V1.size();
  T eSum = 0;
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      eSum += V1(i) * V2(j) * eMat(i, j);
  return eSum;
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

template <typename T> bool IsSymmetricMatrix(MyMatrix<T> const &M) {
  if (M.rows() != M.cols())
    return false;
  int n = M.rows();
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++)
      if (M(i, j) != M(j, i))
        return false;
  return true;
}

template <typename T>
inline typename std::enable_if<!std::is_arithmetic<T>::value, uint32_t>::type
Matrix_Hash(MyMatrix<T> const &M, uint32_t const &seed) {
  std::stringstream s;
  int nbRow = M.rows();
  int nbCol = M.cols();
  for (int iCol = 0; iCol < nbCol; iCol++)
    for (int iRow = 0; iRow < nbRow; iRow++)
      s << " " << M(iRow, iCol);
  std::string converted(s.str());
  const uint8_t *ptr_i = reinterpret_cast<const uint8_t *>(converted.c_str());
  return murmur3_32(ptr_i, converted.size(), seed);
}

template <typename T>
inline typename std::enable_if<std::is_arithmetic<T>::value, uint32_t>::type
Matrix_Hash(MyMatrix<T> const &M, uint32_t const &seed) {
  const T *ptr_T = M.data();
  const uint8_t *ptr_i = reinterpret_cast<const uint8_t *>(ptr_T);
  size_t len = sizeof(T) * M.size();
  return murmur3_32(ptr_i, len, seed);
}

template <typename T>
inline typename std::enable_if<!std::is_arithmetic<T>::value, uint32_t>::type
Vector_Hash(MyVector<T> const &V, uint32_t const &seed) {
  std::stringstream s;
  int n = V.size();
  for (int i = 0; i < n; i++)
    s << " " << V(i);
  std::string converted(s.str());
  const uint8_t *ptr_i = reinterpret_cast<const uint8_t *>(converted.c_str());
  return murmur3_32(ptr_i, converted.size(), seed);
}

template <typename T>
inline typename std::enable_if<std::is_arithmetic<T>::value, uint32_t>::type
Vector_Hash(MyVector<T> const &V, uint32_t const &seed) {
  const T *ptr_T = V.data();
  const uint8_t *ptr_i = reinterpret_cast<const uint8_t *>(ptr_T);
  size_t len = sizeof(T) * V.size();
  return murmur3_32(ptr_i, len, seed);
}

namespace std {
template <typename T> struct hash<MyVector<T>> {
  std::size_t operator()(const MyVector<T> &e_val) const {
    uint32_t seed = 0x1b873540;
    return Vector_Hash(e_val, seed);
  }
};
template <typename T> struct hash<MyMatrix<T>> {
  std::size_t operator()(const MyMatrix<T> &e_val) const {
    uint32_t seed = 0x1b873540;
    return Matrix_Hash(e_val, seed);
  }
};
// clang-format off
}  // namespace std
// clang-format on

template <typename T>
int IntegerDiscriminantInvariant(MyMatrix<T> const &NewMat, int const &n_pes) {
  uint32_t seed = 0x1b873540;
  uint32_t e_hash = Matrix_Hash(NewMat, seed);
  int residue = static_cast<int>(e_hash % n_pes);
  return residue;
}

template <typename T> struct ContainerMatrix {
private:
  MyMatrix<T> const &mat;
  const size_t n_rows;
  const size_t n_cols;
  MyVector<T> v_test;
  std::vector<T> V1, V2;
  std::unordered_set<size_t, std::function<size_t(size_t)>,
                     std::function<bool(size_t, size_t)>>
      set;

public:
  void set_v(std::vector<T> &W, const size_t &idx) {
    if (idx < n_rows) {
      for (size_t i = 0; i < n_cols; i++)
        W[i] = mat(idx, i);
    } else {
      for (size_t i = 0; i < n_cols; i++)
        W[i] = v_test(i);
    }
  }
  ContainerMatrix(MyMatrix<T> const &_mat)
      : mat(_mat), n_rows(mat.rows()), n_cols(mat.cols()), V1(n_cols),
        V2(n_cols) {
    v_test = MyVector<T>(n_cols);
    std::function<size_t(size_t)> fct_hash = [&](size_t idx) -> size_t {
      set_v(V1, idx);
      size_t hash = std::hash<std::vector<T>>()(V1);
      return hash;
    };
    std::function<bool(size_t, size_t)> fct_equal = [&](size_t idx1,
                                                        size_t idx2) -> bool {
      set_v(V1, idx1);
      set_v(V2, idx2);
      for (size_t i = 0; i < n_cols; i++)
        if (V1[i] != V2[i])
          return false;
      return true;
    };
    set = std::unordered_set<size_t, std::function<size_t(size_t)>,
                             std::function<bool(size_t, size_t)>>({}, fct_hash,
                                                                  fct_equal);
    for (size_t i_row = 0; i_row < n_rows; i_row++)
      set.insert(i_row);
  }
  std::optional<size_t> GetIdx() {
    auto iter = set.find(n_rows);
    if (iter == set.end())
      return {};
    size_t idx = *iter;
    return idx;
  }
  std::optional<size_t> GetIdx_v(MyVector<T> const &V) {
    v_test = V;
    return GetIdx();
  }
  template <typename F> std::optional<size_t> GetIdx_f(F f) {
    for (size_t i = 0; i < n_cols; i++)
      v_test(i) = f(i);
    return GetIdx();
  }
};

template <typename T> T GetSmallestMatrixCoefficient(MyMatrix<T> const &M) {
  int nbRow = M.rows();
  int nbCol = M.cols();
  auto get_abs = [](T const &val) -> T {
    if (val < 0)
      return -val;
    return val;
  };
  T the_sma = 1; // Just to shut up warnings. Will not be used
  bool IsAssigned = false;
  auto f_insert = [&](T const &input) -> void {
    T val = get_abs(input);
    if (val > 0) {
      if (!IsAssigned) {
        the_sma = val;
        IsAssigned = true;
      } else {
        if (val < the_sma)
          the_sma = val;
      }
    }
  };
  for (int iCol = 0; iCol < nbCol; iCol++)
    for (int iRow = 0; iRow < nbRow; iRow++)
      f_insert(M(iRow, iCol));
  if (!IsAssigned) {
    std::cerr << "Failed to find a non-zero value for M, so impossible to "
                 "canonicalize\n";
    throw TerminalException{1};
  }
  return the_sma;
}

template <typename T> T GetSmallestVectorCoefficient(MyVector<T> const &V) {
  int n = V.size();
  auto get_abs = [](T const &val) -> T {
    if (val < 0)
      return -val;
    return val;
  };
  T the_sma = 1; // Just to shut up warnings. Will not be used
  bool IsAssigned = false;
  auto f_insert = [&](T const &input) -> void {
    T val = get_abs(input);
    if (val > 0) {
      if (!IsAssigned) {
        the_sma = val;
        IsAssigned = true;
      } else {
        if (val < the_sma)
          the_sma = val;
      }
    }
  };
  for (int i = 0; i < n; i++)
    f_insert(V(i));
  if (!IsAssigned) {
    std::cerr << "Failed to find a non-zero value for V, so impossible to "
                 "canonicalize\n";
    throw TerminalException{1};
  }
  return the_sma;
}

template <typename T> struct ContainerMatrixPositiveScal {
private:
  int nbRow;
  int nbCol;
  std::unique_ptr<ContainerMatrix<T>> Cont;
  MyMatrix<T> Mret;
  MyVector<T> V;

public:
  ContainerMatrixPositiveScal(MyMatrix<T> const &M)
      : nbRow(M.rows()), nbCol(M.cols()) {
    Mret = MyMatrix<T>(nbRow, nbCol);
    V = MyVector<T>(nbCol);
    for (int iRow = 0; iRow < nbRow; iRow++) {
      MyVector<T> V = GetMatrixRow(M, iRow);
      T the_sma = GetSmallestVectorCoefficient(V);
      V /= the_sma;
      AssignMatrixRow(Mret, iRow, V);
    }
    Cont = std::make_unique<ContainerMatrix<T>>(Mret);
  }
  template <typename F> std::optional<size_t> GetIdx_f(F f) {
    for (int i = 0; i < nbCol; i++)
      V(i) = f(i);
    T the_sma = GetSmallestVectorCoefficient(V);
    V /= the_sma;
    return Cont->GetIdx_v(V);
  }
  std::optional<size_t> GetIdx_v(MyVector<T> const &Vin) {
    for (int i = 0; i < nbCol; i++)
      V(i) = Vin(i);
    T the_sma = GetSmallestVectorCoefficient(V);
    V /= the_sma;
    return Cont->GetIdx_v(V);
  }
};

template <typename T>
std::optional<int> get_position_m_v(MyMatrix<T> const &M,
                                    MyVector<T> const &V) {
  int n_rows = M.rows();
  int n_cols = M.cols();
  auto is_in_row = [&](int i_row) -> bool {
    for (int i = 0; i < n_cols; i++) {
      if (M(i_row, i) != V(i)) {
        return false;
      }
    }
    return true;
  };
  for (int i_row = 0; i_row < n_rows; i_row++) {
    if (is_in_row(i_row)) {
      return i_row;
    }
  }
  return {};
}

/*
  G has to be non-degenerate so that we can define the projector.
  --- "G"         is a (n x n) matrix.
  --- "Basis"     is a (p x n) basis of a matrix.
  --- "Basis * G" is a (p x n) matrix representing the (x,u_i) scalar products
  --- "Basis * G * Basis^T" is a matrix of the scalar products
  --- The final matrix is likely to be
      Basis^T (Basis * G * Basis^T) ^ (-1)  Basis G
  If p = n the formula simplifies to Identitity so the formula ought to be
  correct.
 */
template <typename T>
MyMatrix<T> GetProjectionMatrix(MyMatrix<T> const &G,
                                MyMatrix<T> const &Basis) {
  MyMatrix<T> Gred = Basis * G * Basis.transpose();
  if (DeterminantMat(Gred) == 0) {
    std::cerr << "G=\n";
    WriteMatrix(std::cerr, G);
    std::cerr << "Basis=\n";
    WriteMatrix(std::cerr, Basis);
    std::cerr << "Gred=\n";
    WriteMatrix(std::cerr, Gred);
    std::cerr << "The matrix Gred should be invertible\n";
    throw TerminalException{1};
  }
  MyMatrix<T> RetMat = Basis.transpose() * Inverse(Gred) * Basis * G;
  return RetMat;
}

template <typename T>
MyMatrix<T> SubspaceCompletionRational(MyMatrix<T> const &M, int const &n) {
  if (M.rows() == 0)
    return IdentityMat<T>(n);
  return NullspaceTrMat(M);
}

template <typename T> MyVector<T> SignCanonicalizeVector(const MyVector<T> &V) {
  int len = V.size();
  for (int u = 0; u < len; u++) {
    if (V(u) > 0)
      return V;
    if (V(u) < 0)
      return -V;
  }
  std::cerr << "Error in SignCanonicalizeVector\n";
  throw TerminalException{1};
}

template <typename T> MyVector<T> GetDiagonal(MyMatrix<T> const &M) {
  int n_col = M.cols();
  int n_row = M.rows();
  if (n_row != n_col) {
    std::cerr << "The matrix is rectangular but not square, so diagonal is not "
                 "defined\n";
    throw TerminalException{1};
  }
  MyVector<T> V(n_row);
  for (int i = 0; i < n_row; i++) {
    V(i) = M(i, i);
  }
  return V;
}

// clang-format off
#endif  // SRC_MATRIX_MAT_MATRIX_H_
// clang-format on
