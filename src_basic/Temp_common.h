// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_TEMP_COMMON_H_
#define SRC_BASIC_TEMP_COMMON_H_

// Standard includes

// Only standard includes are allowed. Forbidden are:
// ---any boost include
// ---any TBB include
// ---anything that requires a non-trivial installation

// C-style includes

#include <chrono>
#include <ctime>
#include <ctype.h>
#include <getopt.h>
#include <unistd.h>

#include <math.h>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <string>

// Basic C++

#include <exception>

// STL containers

#include <list>
#include <map>
#include <optional>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Functional code

#include <algorithm>
#include <functional>
#include <utility>

// IO streams

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

// Boost serialization

#include "boost_serialization.h"

// Type conversion
#include "Basic_functions.h"
#include "ExceptionsFunc.h"
#include "TypeConversion.h"

// synonyms

typedef unsigned long ulong;
typedef unsigned int uint;

template <typename T> struct is_graphsparseimmutable_class {
  static const bool value = false;
};

template <typename T> struct is_graphbitset_class {
  static const bool value = false;
};

template <typename T> struct is_mymatrix {
  static const bool value = false;
};

unsigned get_random_seed() {
#ifdef USE_NANOSECOND_RAND
  std::timespec ts;
  std::timespec_get(&ts, TIME_UTC);
  unsigned val = ts.tv_nsec;
#else
  unsigned val = time(NULL);
#endif
  return val;
}

void srand_random_set() {
  unsigned val = get_random_seed();
  srand(val);
}

std::string random_string_kernel(std::string const& strChoice, size_t length) {
  const size_t n_index = strChoice.size();
  auto randchar = [&]() -> char {
    size_t pos = size_t(random()) % n_index;
    return strChoice[pos];
  };
  std::string str(length, 0);
  std::generate_n(str.begin(), length, randchar);
  return str;
}

std::string random_string(size_t length) {
  std::string strChoice = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  return random_string_kernel(strChoice, length);
}

std::string random_string_restricted(size_t length) {
  std::string strChoice = "abcdefghijklmnopqrstuvwxyz";
  return random_string_kernel(strChoice, length);
}

std::string GAP_logical(bool const &x) {
  if (x)
    return "true";
  return "false";
}

template <typename T>
void WriteStdVectorStdVectorGAP(std::ostream &os,
                                std::vector<std::vector<T>> const &ListVect) {
  os << "[";
  int IsFirstVect = true;
  for (std::vector<int> const &eVect : ListVect) {
    if (!IsFirstVect)
      os << ",\n";
    IsFirstVect = false;
    os << "[";
    size_t siz = eVect.size();
    for (size_t i = 0; i < siz; i++) {
      if (i > 0)
        os << ",";
      os << eVect[i];
    }
    os << "]";
  }
  os << "]";
}

template <typename T>
void WriteStdVector(std::ostream &os, std::vector<T> const &V) {
  for (auto &eVal : V)
    os << " " << eVal;
  os << "\n";
}

template <typename T>
void WriteStdVectorGAP(std::ostream &os, std::vector<T> const &V) {
  os << "[";
  bool IsFirst = true;
  for (auto &eVal : V) {
    if (!IsFirst)
      os << ",";
    IsFirst = false;
    os << eVal;
  }
  os << "]";
}

template <typename T>
std::istream &operator>>(std::istream &is, std::vector<T> &obj) {
  int n;
  is >> n;
  obj.resize(n);
  for (int i = 0; i < n; i++) {
    T eVal;
    is >> eVal;
    obj[i] = eVal;
  }
  return is;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, std::vector<T> const &ListVal) {
  int n = ListVal.size();
  os << " " << n;
  for (int i = 0; i < n; i++)
    os << " " << ListVal[i];
  return os;
}

template <typename T> struct CollectedResult {
  std::vector<T> LVal;
  std::vector<int> LMult;
};

template <typename T>
CollectedResult<T> Collected(std::vector<T> const &eVect) {
  std::set<T> SetVal;
  for (auto &eVal : eVect)
    SetVal.insert(eVal);
  std::vector<T> LVal;
  for (auto &eVal : SetVal)
    LVal.push_back(eVal);
  size_t eSize = LVal.size();
  std::vector<int> LMult(eSize, 0);
  auto UpPosition = [&](T const &eVal) -> void {
    for (size_t i = 0; i < eSize; i++)
      if (LVal[i] == eVal) {
        LMult[i] += 1;
        return;
      }
    std::cerr << "Should never reach that stage\n";
    throw TerminalException{1};
  };
  for (auto &eVal : eVect)
    UpPosition(eVal);
  return {std::move(LVal), std::move(LMult)};
}

std::vector<int> StdVectorFirstNentries(size_t const &N) {
  std::vector<int> eList(N);
  for (size_t i = 0; i < N; i++)
    eList[i] = static_cast<int>(i);
  return eList;
}

template <typename T>
void WriteVectorInt_GAP(std::ostream &os, std::vector<T> const &OneInc) {
  size_t siz = OneInc.size();
  os << "[";
  for (size_t i = 0; i < siz; i++) {
    if (i > 0)
      os << ",";
    size_t eVal = OneInc[i] + 1;
    os << eVal;
  }
  os << "]";
}

std::vector<int> DivideListPosition(int const &len, int const &nbBlock) {
  std::vector<int> ListVal;
  double fact = static_cast<double>(len) / static_cast<double>(nbBlock);
  for (int i = 0; i <= nbBlock; i++) {
    double pos_d = fact * static_cast<double>(i);
    int pos_i;
    NearestInteger_double_int(pos_d, pos_i);
    pos_i = std::max(0, pos_i);
    pos_i = std::min(len, pos_i);
    ListVal.push_back(pos_i);
  }
  return ListVal;
}

// Returned fields {v1,v2}
// v1 is the field returning the sorted index to the original index
// v2 is the field returning the original index to the sorted index
template <typename T>
std::pair<std::vector<int>, std::vector<int>>
SortingLists(std::vector<T> const &ListV) {
  struct PairData {
    size_t i;
    T x;
  };
  std::size_t len = ListV.size();
  std::vector<PairData> ListPair(len);
  for (std::size_t i = 0; i < len; i++) {
    PairData ePair{i, ListV[i]};
    ListPair[i] = ePair;
  }
  sort(ListPair.begin(), ListPair.end(),
       [](PairData const &x1, PairData const &x2) -> bool {
         if (x1.x < x2.x)
           return true;
         if (x2.x < x1.x)
           return false;
         return x1.i < x2.i;
       });
  std::vector<int> v1(len);
  std::vector<int> v2(len);
  for (size_t i = 0; i < len; i++) {
    size_t eIdx = static_cast<int>(ListPair[i].i);
    v1[i] = static_cast<int>(eIdx);
    v2[eIdx] = static_cast<int>(i);
  }
  return {std::move(v1), std::move(v2)};
}

inline int T_NormGen(int const &x) { return abs(x); }

inline int CanonicalizationUnit(int const &eVal) {
  if (eVal < 0)
    return -1;
  return 1;
}

// clang-format off
#endif  // SRC_BASIC_TEMP_COMMON_H_
// clang-format on
