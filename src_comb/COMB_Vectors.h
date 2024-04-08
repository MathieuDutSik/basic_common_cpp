// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_COMB_COMB_VECTORS_H_
#define SRC_COMB_COMB_VECTORS_H_

#include <set>
#include <unordered_set>
#include <vector>

template <typename T>
std::vector<T> UnionVect(std::vector<T> const &V1, std::vector<T> const &V2) {
  std::unordered_set<T> eSet;
  for (auto &eVal : V1)
    eSet.insert(eVal);
  for (auto &eVal : V2)
    eSet.insert(eVal);
  std::vector<T> eV;
  for (auto &eVal : eSet)
    eV.push_back(eVal);
  return eV;
}
template <typename T>
std::vector<T> DiffIntVect(std::vector<T> const &V1, std::vector<T> const &V2,
                           int n_int) {
  std::unordered_set<T> eSet2;
  for (auto &eVal : V2)
    eSet2.insert(eVal);
  std::vector<T> eV;
  for (auto &eVal : V1) {
    if (eSet2.count(eVal) == static_cast<size_t>(n_int))
      eV.push_back(eVal);
  }
  return eV;
}

// The difference V1 - V2
template <typename T>
std::vector<T> DifferenceVect(std::vector<T> const &V1,
                              std::vector<T> const &V2) {
  return DiffIntVect(V1, V2, 0);
}

// The intersection V1 \cap V2
template <typename T>
std::vector<T> IntersectionVect(std::vector<T> const &V1,
                                std::vector<T> const &V2) {
  return DiffIntVect(V1, V2, 1);
}

template <typename T>
void AppendVect(std::vector<T> &V1, std::vector<T> const &V2) {
  for (auto &eVal : V2)
    V1.push_back(eVal);
}

template <typename T>
bool IsSubset(std::vector<T> const &S1, std::vector<T> const &S2) {
  for (auto &eVal : S2) {
    int pos = PositionVect(S1, eVal);
    if (pos == -1)
      return false;
  }
  return true;
}

// V1 and V2 should have no repetition
template<typename T>
bool IsEqualSet(std::vector<T> const& V1, std::vector<T> const& V2) {
  if (V1.size() != V2.size()) {
    return false;
  }
  std::unordered_set<T> set;
  for (auto &val : V1) {
    set.insert(val);
  }
  for (auto & val : V2) {
    if (set.count(val) == 0) {
      return false;
    }
  }
  return true;
}

template <typename T> std::vector<T> VectorAsSet(std::vector<T> const &V) {
  std::unordered_set<T> eSet;
  for (auto &eVal : V)
    eSet.insert(eVal);
  std::vector<T> eV;
  for (auto &eVal : eSet)
    eV.push_back(eVal);
  return eV;
}

template <typename T> struct popable_vector {
private:
  std::vector<T> V;
  size_t pos;

public:
  popable_vector() : pos(0) {}
  popable_vector(std::vector<T> const &_V) : V(_V), pos(_V.size()) {}
  size_t size() const { return pos; }
  void push_back(T const &val) {
    if (V.size() == pos) {
      V.push_back(val);
      pos++;
      return;
    }
    V[pos] = val;
    pos++;
  }
  T pop() {
    pos--;
    return V[pos];
  }
};

// clang-format off
#endif  // SRC_COMB_COMB_VECTORS_H_
// clang-format on
