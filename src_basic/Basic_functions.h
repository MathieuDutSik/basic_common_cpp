// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_BASIC_BASIC_FUNCTIONS_H_
#define SRC_BASIC_BASIC_FUNCTIONS_H_

#include <algorithm>
#include <iostream>
#include <map>
#include <ranges>
#include <sstream>
#include <string>
#include <vector>

template <typename T> T VectorSum(std::vector<T> const &eVect) {
  T eSum = 0;
  for (T eVal : eVect)
    eSum += eVal;
  return eSum;
}

template <typename T> T VectorMin(std::vector<T> const &eVect) {
  T eMin = eVect[0];
  for (T eVal : eVect)
    if (eVal < eMin)
      eMin = eVal;
  return eMin;
}

template <typename T> T VectorMax(std::vector<T> const &eVect) {
  T eMax = eVect[0];
  for (T eVal : eVect)
    if (eVal > eMax)
      eMax = eVal;
  return eMax;
}

template <typename T> T VectorAvg(std::vector<T> const &eVect) {
  T eSum = 0;
  for (T eVal : eVect)
    eSum += eVal;
  T eAvg = eSum / static_cast<int>(eVect.size());
  return eAvg;
}

template <typename T>
void ShowAttainmentVector(std::ostream &os, std::vector<T> const &eVect) {
  std::map<T, int> TheMap;
  for (auto &eVal : eVect) {
    auto search = TheMap.find(eVal);
    if (search == TheMap.end()) {
      TheMap[eVal] = 1;
    } else {
      TheMap[eVal]++;
    }
  }
  for (auto &ePair : TheMap)
    os << "Value " << ePair.first << " attained " << ePair.second << " times\n";
}

// use std::max for generic types (int, long, float, ...)
template <typename T> constexpr T T_max(T const &eVal1, T const &eVal2) {
  if (eVal1 > eVal2)
    return eVal1;
  return eVal2;
}

template <typename T> constexpr T T_abs(T const &eVal) {
  if (eVal > 0)
    return eVal;
  return -eVal;
}

// use std::min for generic types (int, long, float, ...)
template <typename T> constexpr T T_min(T const &eVal1, T const &eVal2) {
  if (eVal1 > eVal2)
    return eVal2;
  return eVal1;
}

template <typename T> constexpr int T_sign(T const &eVal) {
  if (eVal > 0)
    return 1;
  if (eVal < 0)
    return -1;
  return 0;
}

template <typename T> constexpr T NextIdx(T const &len, T const &i) {
  if (i == len - 1)
    return 0;
  return i + 1;
}

template <typename T> constexpr T PrevIdx(T const &len, T const &i) {
  if (i == 0)
    return len - 1;
  return i - 1;
}

template <typename T> constexpr T MyPow(T const &eVal, int const &n) {
  T eRet = 1;
  for (int i = 0; i < n; i++)
    eRet *= eVal;
  return eRet;
}

template <typename T> int PositionVect(std::vector<T> const &V, T const &eVal) {
  auto it = std::ranges::find(V, eVal);
  if (it == V.end())
    return -1;
  return static_cast<int>(it - V.begin());
}

template <typename T> bool IsVectorConstant(std::vector<T> const &V) {
  if (V.empty())
    return true;
  return std::ranges::all_of(V, [&](auto const &x) { return x == V[0]; });
}

template <typename T>
std::vector<T> ConcatenationTwo(std::vector<T> const &L1,
                                std::vector<T> const &L2) {
  std::vector<T> ret = L1;
  ret.insert(ret.end(), L2.begin(), L2.end());
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

template <typename T> std::vector<T> Concatenation(std::vector<T> const &L) {
  return L;
}

template <typename T, typename... Args>
std::vector<T> Concatenation(std::vector<T> const &first, Args... args) {
  std::vector<T> ret = first;
  for (auto &eVal : Concatenation(args...))
    ret.push_back(eVal);
  return ret;
}

template <typename T, class UnaryPredicate>
std::vector<T> ListT(std::vector<T> const &V, UnaryPredicate const &f) {
  std::vector<T> retV;
  for (auto &eVal : V)
    retV.push_back(f(eVal));
  return retV;
}

template <typename T, class UnaryPredicate>
int PositionProperty(std::vector<T> const &V, UnaryPredicate const &f) {
  int len = V.size();
  for (int i = 0; i < len; i++)
    if (f(V[i]))
      return i;
  return -1;
}

template <typename T, class UnaryPredicate>
bool ForAll(std::vector<T> const &V, UnaryPredicate const &f) {
  return std::ranges::all_of(V, f);
}

template <typename T, class UnaryPredicate>
std::vector<T> Filtered(std::vector<T> const &V, UnaryPredicate const &f) {
  std::vector<T> LRet;
  for (auto &eVal : V)
    if (f(eVal))
      LRet.push_back(eVal);
  return LRet;
}

// Parsing strings (It takes a std::string because we would need a const char*
// with a null terminated string and this we would not have with string_view

template <typename T> T ParseScalar(std::string const &estr) {
  T ret_val;
  std::istringstream(estr) >> ret_val;
  return ret_val;
}

template <typename T>
void ParseScalar_inplace(std::string const &estr, T &ret_val) {
  std::istringstream(estr) >> ret_val;
}

// clang-format off
#endif  // SRC_BASIC_BASIC_FUNCTIONS_H_
// clang-format on
