// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_COMB_BOOST_BITSET_H_
#define SRC_COMB_BOOST_BITSET_H_

// Boost libraries

#include "Boost_bitset_kernel.h"
#include "Temp_common.h"
#include <string>
#include <vector>

std::vector<int> FaceTo01vector(Face const &eSet) {
  size_t nbVert = eSet.size();
  size_t siz = eSet.count();
  std::vector<int> eList(nbVert, 0);
  boost::dynamic_bitset<>::size_type aRow = eSet.find_first();
  for (size_t i = 0; i < siz; i++) {
    eList[aRow] = 1;
    aRow = eSet.find_next(aRow);
  }
  return eList;
}

template <typename Tidx>
Face VectorToFace(std::vector<Tidx> const &V, size_t siz) {
  Face f(siz);
  for (auto &eVal : V) {
    f[eVal] = 1;
  }
  return f;
}

std::string StringFace(Face const &f) {
  std::string str_ret;
  size_t len = f.size();
  for (size_t i = 0; i < len; i++)
    str_ret += std::to_string(static_cast<int>(f[i]));
  return str_ret;
}

std::string SignatureFace(Face const &f) {
  std::string str_ret =
      std::to_string(f.count()) + "/" + std::to_string(f.size());
  return str_ret;
}

void WriteFace(std::ostream &os, Face const &eList) {
  size_t len = eList.size();
  os << len;
  for (size_t i = 0; i < len; i++) {
    int eVal = eList[i];
    os << " " << eVal;
  }
  os << "\n";
}

Face ReadFace(std::istream &is) {
  if (!is.good()) {
    std::cerr << "ReadFace operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  size_t len;
  int eVal;
  is >> len;
  Face eFace(len);
  for (size_t i = 0; i < len; i++) {
    is >> eVal;
    eFace[i] = eVal;
  }
  return eFace;
}

Face ReadFaceFile(std::string const &eFile) {
  std::ifstream is(eFile);
  return ReadFace(is);
}

Face RandomFace(int n) {
  Face f(n);
  for (int i = 0; i < n; i++) {
    int eVal = random() % 2;
    f[i] = eVal;
  }
  return f;
}

Face RandomKFace(int n, int k) {
  if (n < 0) {
    std::cerr << "We should have n >= 0. We have n=" << n << "\n";
    throw TerminalException{1};
  }
  if (k < 0 || k > n) {
    std::cerr << "We should have 0 <= k <= n. We have n=" << n << " k=" << k
              << "\n";
    throw TerminalException{1};
  }
  Face f(n);
  if (2 * k < n) {
    int n_done = 0;
    while (true) {
      int pos = random() % n;
      if (f[pos] == 0) {
        f[pos] = 1;
        n_done++;
      }
      if (n_done == k)
        break;
    }
  } else {
    for (int i = 0; i < n; i++)
      f[i] = 0;
    int n_done = n;
    while (true) {
      int pos = random() % n;
      if (f[pos] == 1) {
        f[pos] = 0;
        n_done--;
      }
      if (n_done == k)
        break;
    }
  }
  return f;
}

vectface ReadListFace(std::istream &is) {
  if (!is.good()) {
    std::cerr << "ReadListFace operation failed because stream is not valid\n";
    throw TerminalException{1};
  }
  size_t nbFace;
  is >> nbFace;
  if (nbFace == 0) {
    std::cerr << "We cannot handle that case because we need the base length\n";
    throw TerminalException{1};
  }
  Face f = ReadFace(is);
  vectface ListFace(f.size());
  ListFace.push_back(f);
  for (size_t iFace = 1; iFace < nbFace; iFace++) {
    Face f2 = ReadFace(is);
    ListFace.push_back(f2);
  }
  return ListFace;
}

void WriteListFace(std::ostream &os, vectface const &ListFace) {
  size_t nbFace = ListFace.size();
  os << nbFace << "\n";
  for (size_t iFace = 0; iFace < nbFace; iFace++)
    WriteFace(os, ListFace[iFace]);
}

void WriteFaceBracket(std::ostream &os, Face const &f, int shift) {
  size_t nb = f.count();
  //  int siz=f.size();
  os << "[";
  boost::dynamic_bitset<>::size_type aPos = f.find_first();
  for (size_t i = 0; i < nb; i++) {
    if (i > 0)
      os << ",";
    int eVal = static_cast<int>(aPos) + shift;
    os << eVal;
    aPos = f.find_next(aPos);
  }
  os << "]";
}

void WriteFaceGAP(std::ostream &os, Face const &f) {
  int shift = 1;
  WriteFaceBracket(os, f, shift);
}

void WriteFacePYTHON(std::ostream &os, Face const &f) {
  int shift = 0;
  WriteFaceBracket(os, f, shift);
}

void WriteListFaceGAP(std::ostream &os, vectface const &ListFace) {
  os << "[";
  bool IsFirst = true;
  for (auto &eFace : ListFace) {
    if (!IsFirst)
      os << ",";
    IsFirst = false;
    WriteFaceGAP(os, eFace);
  }
  os << "]";
}

void WriteListFacePYTHON(std::ostream &os, vectface const &ListFace) {
  os << "[";
  bool IsFirst = true;
  for (auto &eFace : ListFace) {
    if (!IsFirst)
      os << ",";
    IsFirst = false;
    WriteFacePYTHON(os, eFace);
  }
  os << "]";
}

void WriteListFaceGAPfile(std::string const &eFile, vectface const &ListFace) {
  std::ofstream os(eFile);
  os << "return ";
  WriteListFaceGAP(os, ListFace);
  os << ";\n";
}

// We require x and y to be of the same size
bool operator<(Face const &x, Face const &y) {
  size_t len = x.size();
  for (size_t i = 0; i < len; i++) {
    if (x[i] == 0 && y[i] == 1)
      return true;
    if (x[i] == 1 && y[i] == 0)
      return false;
  }
  return false;
}

void PrintVectInt(std::ostream &os, Face const &eList) {
  size_t len = eList.size();
  for (size_t i = 0; i < len; i++)
    if (eList[i] == 1)
      os << " " << i;
  os << "\n";
}

Face FullFace(size_t const &len) {
  Face eFace(len);
  for (size_t u = 0; u < len; u++)
    eFace[u] = 1;
  return eFace;
}

ulong FaceToUnsignedLong(Face const &f) {
  size_t len = f.size();
  if (len > 32) {
    std::cerr << "Too large value, conversion impossible";
    throw TerminalException{1};
  }
  ulong pos = 0;
  ulong pow = 1;
  for (size_t i = 0; i < len; i++) {
    pos += pow * f[i];
    pow *= 2;
  }
  return pos;
}

Face UnsignedLongToFace(size_t const &len, ulong const &eVal) {
  if (len > 32) {
    std::cerr << "length error\n";
    throw TerminalException{1};
  }
  ulong eWork = eVal;
  Face eFace(len);
  ulong pow = 1;
  for (size_t i = 0; i < len; i++) {
    ulong Pow2 = pow * 2;
    ulong res = eWork % Pow2;
    if (res == pow) {
      eFace[i] = 1;
      eWork -= pow;
    }
    pow = Pow2;
  }
  return eFace;
}

void VectVectInt_Print_Kernel(std::ostream &os, vectface const &ListOrbit,
                              int const &shift) {
  size_t nbOrbit = ListOrbit.size();
  os << "[";
  for (size_t iOrbit = 0; iOrbit < nbOrbit; iOrbit++) {
    if (iOrbit > 0)
      os << ",\n";
    Face eRepr = ListOrbit[iOrbit];
    size_t siz = eRepr.count();
    os << "[";
    boost::dynamic_bitset<>::size_type eVal = eRepr.find_first();
    for (size_t i = 0; i < siz; i++) {
      if (i > 0)
        os << ",";
      os << int(eVal + shift);
      eVal = eRepr.find_next(eVal);
    }
    os << "]";
  }
  os << "]";
}

void VectVectInt_Magma_Print(std::ostream &os, vectface const &ListOrbit) {
  VectVectInt_Print_Kernel(os, ListOrbit, 0);
}

void VectVectInt_Magma_PrintFile(std::string const &eFile,
                                 vectface const &ListOrbit) {
  std::ofstream os(eFile);
  os << "return ";
  VectVectInt_Magma_Print(os, ListOrbit);
  os << ";\n";
}

void VectVectInt_Gap_Print(std::ostream &os, vectface const &ListOrbit) {
  VectVectInt_Print_Kernel(os, ListOrbit, 1);
}

void VectVectInt_Gap_PrintFile(std::string const &eFile,
                               vectface const &ListOrbit) {
  std::ofstream os(eFile);
  os << "return ";
  VectVectInt_Gap_Print(os, ListOrbit);
  os << ";\n";
}

template <typename T>
void VectVectInt_SetInt_PrintFile(std::string const &eFile,
                                  vectface const &ListOrbit) {
  std::ofstream os(eFile);
  os << ListOrbit.size() << "\n";
  for (const Face &f : ListOrbit) {
    T res = getsetasint<T>(f);
    os << res << "\n";
  }
}

// clang-format off
#endif  // SRC_COMB_BOOST_BITSET_H_
// clang-format on
