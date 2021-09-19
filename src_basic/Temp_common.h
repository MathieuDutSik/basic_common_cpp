#ifndef TEMP_COMMON_TO_ALL
#define TEMP_COMMON_TO_ALL

// Standard includes

// All the code here should have common usage with
// ---Oceanography
// ---Chemistry
// ---Mathematics
// If it applies to only one of those fields, then it should
// be out.

// Only standard includes are allowed. Forbidden are:
// ---any boost include
// ---any TBB include
// ---anything that requires a non-trivial installation

// C-style includes

#include <ctype.h>
#include <unistd.h>
#include <getopt.h>
#include <chrono>
#include <ctime>


#include <math.h>

#include <string>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

// Basic C++

#include <exception>

// STL containers

#include <vector>
#include <list>
#include <set>
#include <unordered_set>
#include <map>
#include <unordered_map>

// Functional code

#include <functional>
#include <algorithm>

// IO streams

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

// Boost serialization

#include <boost/archive/tmpdir.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/assume_abstract.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/split_free.hpp>


// Type conversion
#include "TypeConversion.h"
#include "ExceptionEnding.h"

// synonyms

typedef unsigned long ulong;
typedef unsigned int uint;


template<typename T>
struct is_graphsparseimmutable_class {
  static const bool value = false;
};

template<typename T>
struct is_graphbitset_class {
  static const bool value = false;
};



template<typename T>
T VectorSum(std::vector<T> const& eVect)
{
  T eSum=0;
  for (T eVal : eVect)
    eSum += eVal;
  return eSum;
}

template<typename T>
T VectorMin(std::vector<T> const& eVect)
{
  T eMin=eVect[0];
  for (T eVal : eVect)
    if (eVal < eMin)
      eMin=eVal;
  return eMin;
}

template<typename T>
T VectorMax(std::vector<T> const& eVect)
{
  T eMax=eVect[0];
  for (T eVal : eVect)
    if (eVal > eMax)
      eMax=eVal;
  return eMax;
}

template<typename T>
T VectorAvg(std::vector<T> const& eVect)
{
  T eSum=0;
  for (T eVal : eVect)
    eSum += eVal;
  T eAvg = eSum / int(eVect.size());
  return eAvg;
}




template<typename T>
void ShowAttainmentVector(std::ostream & os, std::vector<T> const& eVect)
{
  std::map<T,int> TheMap;
  for (auto & eVal : eVect) {
    auto search = TheMap.find(eVal);
    if (search == TheMap.end()) {
      TheMap[eVal] = 1;
    }
    else {
      TheMap[eVal]++;
    }
  }
  for (auto & ePair : TheMap)
    os << "Value " << ePair.first << " attained " << ePair.second << " times\n";
}



// use std::max for generic types (int, long, float, ...)
template<typename T>
T T_max(T const& eVal1, T const& eVal2)
{
  if (eVal1 > eVal2)
    return eVal1;
  return eVal2;
}


// use std::min for generic types (int, long, float, ...)
template<typename T>
T T_min(T const& eVal1, T const& eVal2)
{
  if (eVal1 > eVal2)
    return eVal2;
  return eVal1;
}


template<typename T>
int T_sign(T const& eVal)
{
  if (eVal > 0)
    return 1;
  if (eVal < 0)
    return -1;
  return 0;
}

void srand_random_set()
{
#ifdef USE_NANOSECOND_RAND
  std::timespec ts;
  std::timespec_get(&ts, TIME_UTC);
  unsigned val = ts.tv_nsec;
  srand(val);
#else
  srand(time(NULL));
#endif
}



std::string random_string( size_t length )
{
  srand_random_set();
  auto randchar = []() -> char {
    const char charset[] =
    "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";
    const size_t max_index = (sizeof(charset) - 1);
    return charset[ size_t(rand()) % max_index ];
  };
  std::string str(length,0);
  std::generate_n( str.begin(), length, randchar );
  return str;
}



std::string random_string_restricted( size_t length )
{
  srand ( time(NULL) );
  auto randchar = []() -> char {
    const char charset[] = "abcdefghijklmnopqrstuvwxyz";
    const size_t max_index = (sizeof(charset) - 1);
    return charset[ size_t(rand()) % max_index ];
  };
  std::string str(length,0);
  std::generate_n( str.begin(), length, randchar );
  return str;
}


std::string GAP_logical(bool const& x)
{
  if (x)
    return "true";
  return "false";
}



template<typename T>
void WriteStdVectorStdVectorGAP(std::ostream & os, std::vector<std::vector<T> > const& ListVect)
{
  os << "[";
  int IsFirstVect=true;
  for (std::vector<int> const& eVect : ListVect) {
    if (!IsFirstVect)
      os << ",\n";
    IsFirstVect=false;
    os << "[";
    size_t siz=eVect.size();
    for (size_t i=0; i<siz; i++) {
      if (i>0)
        os << ",";
      os << eVect[i];
    }
    os << "]";
  }
  os << "]";
}

template<typename T>
void WriteStdVector(std::ostream& os, std::vector<T> const& V)
{
  for (auto & eVal : V)
    os << " " << eVal;
  os << "\n";
}


template<typename T>
void WriteStdVectorGAP(std::ostream& os, std::vector<T> const& V)
{
  os << "[";
  bool IsFirst=true;
  for (auto & eVal : V) {
    if (!IsFirst)
      os << ",";
    IsFirst=false;
    os << eVal;
  }
  os << "]";
}

template<typename T>
std::istream& operator>>(std::istream& is, std::vector<T>& obj)
{
  int n;
  is >> n;
  obj.resize(n);
  for (int i=0; i<n; i++) {
    T eVal;
    is >> eVal;
    obj[i]=eVal;
  }
  return is;
}
template<typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T> const& ListVal)
{
  int n=ListVal.size();
  os << " " << n;
  for (int i=0; i<n; i++)
    os << " " << ListVal[i];
  return os;
}






template<typename T>
struct CollectedResult {
  std::vector<T> LVal;
  std::vector<int> LMult;
};

template<typename T>
CollectedResult<T> Collected(std::vector<T> const& eVect)
{
  std::set<T> SetVal;
  for (auto & eVal : eVect)
    SetVal.insert(eVal);
  std::vector<T> LVal;
  for (auto & eVal : SetVal)
    LVal.push_back(eVal);
  size_t eSize=LVal.size();
  std::vector<int> LMult(eSize,0);
  auto UpPosition=[&](T const& eVal) -> void {
    for (size_t i=0; i<eSize; i++)
      if (LVal[i] == eVal) {
	LMult[i] += 1;
	return;
      }
    std::cerr << "Should never reach that stage\n";
    throw TerminalException{1};
  };
  for (auto & eVal : eVect)
    UpPosition(eVal);
  return {std::move(LVal), std::move(LMult)};
}

template<typename T>
inline T NextIdx(T const& len, T const& i)
{
  if (i == len-1)
    return 0;
  return i+1;
}

template<typename T>
inline T PrevIdx(T const& len, T const& i)
{
  if (i == 0)
    return len-1;
  return i-1;
}


template<typename T>
T MyPow(T const& eVal, int const& n)
{
  T eRet=1;
  for (int i=0; i<n; i++)
    eRet *= eVal;
  return eRet;
}


std::vector<int> StdVectorFirstNentries(size_t const& N)
{
  std::vector<int> eList(N);
  for (size_t i=0; i<N; i++)
    eList[i] = int(i);
  return eList;
}

template<typename T>
int PositionVect(std::vector<T> const& V, T const& eVal)
{
  size_t len=V.size();
  for (size_t i=0; i<len; i++)
    if (V[i] == eVal)
      return int(i);
  return -1;
}

template<typename T>
bool IsVectorConstant(std::vector<T> const& V)
{
  if (V.size() == 0)
    return true;
  T eVal=V[0];
  for (auto & fVal : V) {
    if (eVal != fVal)
      return false;
  }
  return true;
}


template<typename T>
void WriteVectorInt_GAP(std::ostream &os, std::vector<T> const& OneInc)
{
  size_t siz=OneInc.size();
  os << "[";
  for (size_t i=0; i<siz; i++) {
    if (i>0)
      os << ",";
    size_t eVal=OneInc[i]+1;
    os << eVal;
  }
  os << "]";
}


std::vector<int> DivideListPosition(int const& len, int const& nbBlock)
{
  std::vector<int> ListVal;
  for (int i=0; i<=nbBlock; i++) {
    double pos_d = (double(i) / double(nbBlock)) * double(len);
    int pos_i;
    NearestInteger_double_int(pos_d, pos_i);
    pos_i = std::max(0, pos_i);
    pos_i = std::min(len, pos_i);
    ListVal.push_back(pos_i);
  }
  return ListVal;
}


template<typename T>
std::vector<T> ConcatenationTwo(std::vector<T> const& L1, std::vector<T> const& L2)
{
  std::vector<T> ret = L1;
  for (auto & eVal : L2)
    ret.push_back(eVal);
  return ret;
}



/*
template<typename U>
std::vector<U> operator+(std::vector<U> const& L1, std::vector<U> const& L2)
{
  std::vector<U> ret = L1;
  for (auto & eVal : L2)
    ret.push_back(eVal);
  return ret;
}
*/



template<typename T>
std::vector<T> Concatenation(std::vector<T> const& L)
{
  return L;
}


template<typename T, typename... Args>
std::vector<T> Concatenation(std::vector<T> const& first, Args... args)
{
  std::vector<T> ret = first;
  for (auto & eVal : Concatenation(args...))
    ret.push_back(eVal);
  return ret;
}


template<typename T, class UnaryPredicate>
std::vector<T> ListT(std::vector<T> const& V, UnaryPredicate const& f)
{
  std::vector<T> retV;
  for (auto & eVal : V)
    retV.push_back(f(eVal));
  return retV;
}




template<typename T, class UnaryPredicate>
int PositionProperty(std::vector<T> const& V, UnaryPredicate const& f)
{
  int len=V.size();
  for (int i=0; i<len; i++)
    if (f(V[i]))
      return i;
  return -1;
}

// Actually there is a version all_of in <algorithm> but it uses iterators
// which is kind of painful.

template<typename T, class UnaryPredicate>
bool ForAll(std::vector<T> const& V, UnaryPredicate const& f)
{
  for (auto & eVal : V)
    if (!f(eVal))
      return false;
  return true;
}

template<typename T, class UnaryPredicate>
std::vector<T> Filtered(std::vector<T> const& V, UnaryPredicate const& f)
{
  std::vector<T> LRet;
  for (auto & eVal : V)
    if (f(eVal))
      LRet.push_back(eVal);
  return LRet;
}

// Returned fields {v1,v2}
// v1 is the field returning the sorted index to the original index
// v2 is the field returning the original index to the sorted index
template<typename T>
std::pair<std::vector<int>,std::vector<int>> SortingLists(std::vector<T> const & ListV)
{
  struct PairData {
    size_t i;
    T x;
  };
  std::size_t len=ListV.size();
  std::vector<PairData> ListPair(len);
  for (std::size_t i=0; i<len; i++) {
    PairData ePair{i, ListV[i]};
    ListPair[i]=ePair;
  }
  sort(ListPair.begin(), ListPair.end(),
       [](PairData const & x1, PairData const& x2) -> bool {
         if (x1.x < x2.x)
           return true;
         if (x2.x < x1.x)
           return false;
         return x1.i< x2.i;
       });
  std::vector<int> v1(len);
  std::vector<int> v2(len);
  for (size_t i=0; i<len; i++) {
    size_t eIdx=int(ListPair[i].i);
    v1[i]=int(eIdx);
    v2[eIdx]=int(i);
  }
  return {std::move(v1), std::move(v2)};
}



#endif
