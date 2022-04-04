#ifndef SRC_BASIC_BASIC_FUNCTIONS_H_
#define SRC_BASIC_BASIC_FUNCTIONS_H_

#include <vector>
#include <map>
#include <iostream>


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


template<typename T>
T T_abs(T const& eVal)
{
  if (eVal > 0)
    return eVal;
  T fVal= - eVal;
  return fVal;
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

// Parsing strings (It takes a std::string because we would need a const char* with a null terminated string
// and this we would not have with string_view

template<typename T>
T ParseScalar(std::string const& estr)
{
  T ret_val;
  std::istringstream is(estr);
  is >> ret_val;
  return ret_val;
}

template<typename T>
void ParseScalar_inplace(std::string const& estr, T & ret_val)
{
  std::istringstream is(estr);
  is >> ret_val;
}


#endif // SRC_BASIC_BASIC_FUNCTIONS_H_
