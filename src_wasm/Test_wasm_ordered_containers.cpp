// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// Wasm test: MyVector / MyMatrix used as keys of the ordered containers and
// as elements of the ordering algorithms.
//
// This exists because of a failure that no other test could see. MyVector<T>
// is an alias of Eigen::Matrix, a type this project does not define, so the
// only lookup that reaches an operator of it from inside the standard library
// is the argument dependent one. Argument dependent lookup searches the
// namespace of Eigen and the namespaces of the template arguments, so an
// operator< at global scope is found when T is at global scope (mpq_class of
// gmpxx) and is not found when T is in a namespace
// (boost::multiprecision::cpp_rational). The specialization of std::less does
// not save it: libc++ 220108 rewrites the comparator std::less<Key> of an
// ordered container into its transparent form and compares with a < b, which
// [comparisons.less] entitles it to do, so the specialization is never
// instantiated.
//
// The scalar types below are therefore chosen on purpose: cpp_rational and
// cpp_int live in a namespace, which is the case that fails.

#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryCommon.h"
#include "MAT_Matrix.h"
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>

using cpp_int = boost::multiprecision::cpp_int;
using cpp_rational = boost::multiprecision::cpp_rational;

static int n_fail = 0;

#define CHECK(cond)                                                            \
  do {                                                                         \
    if (!(cond)) {                                                             \
      std::cerr << "FAIL " << __FILE__ << ":" << __LINE__ << " : " << #cond    \
                << "\n";                                                       \
      n_fail++;                                                                \
    }                                                                          \
  } while (0)

template <typename T> MyVector<T> get_vector(int a, int b) {
  MyVector<T> V(2);
  V(0) = a;
  V(1) = b;
  return V;
}

template <typename T> MyMatrix<T> get_matrix(int a, int b, int c, int d) {
  MyMatrix<T> M(2, 2);
  M(0, 0) = a;
  M(0, 1) = b;
  M(1, 0) = c;
  M(1, 1) = d;
  return M;
}

// std::set of vectors: needs the ordering from inside the standard library.
template <typename T> void test_set_vector() {
  std::set<MyVector<T>> S;
  S.insert(get_vector<T>(1, 2));
  S.insert(get_vector<T>(1, 2));
  S.insert(get_vector<T>(0, 5));
  CHECK(S.size() == 2);
  CHECK(S.count(get_vector<T>(1, 2)) == 1);
  CHECK(S.count(get_vector<T>(7, 7)) == 0);
  // The smallest is the lexicographic one.
  CHECK(*S.begin() == get_vector<T>(0, 5));
}

// std::map keyed by a vector.
template <typename T> void test_map_vector() {
  std::map<MyVector<T>, int> M;
  M[get_vector<T>(1, 2)] = 3;
  M[get_vector<T>(1, 2)] = 4;
  M[get_vector<T>(2, 1)] = 5;
  CHECK(M.size() == 2);
  CHECK(M[get_vector<T>(1, 2)] == 4);
}

// std::set of matrices.
template <typename T> void test_set_matrix() {
  std::set<MyMatrix<T>> S;
  S.insert(get_matrix<T>(1, 0, 0, 1));
  S.insert(get_matrix<T>(1, 0, 0, 1));
  S.insert(get_matrix<T>(0, 1, 1, 0));
  CHECK(S.size() == 2);
}

// std::sort never consults std::less<Key>: it compares with a < b directly.
// This is the case the container aliases would not have covered.
template <typename T> void test_sort_vector() {
  std::vector<MyVector<T>> L;
  L.push_back(get_vector<T>(2, 0));
  L.push_back(get_vector<T>(0, 1));
  L.push_back(get_vector<T>(1, 9));
  std::sort(L.begin(), L.end());
  CHECK(L[0] == get_vector<T>(0, 1));
  CHECK(L[1] == get_vector<T>(1, 9));
  CHECK(L[2] == get_vector<T>(2, 0));
  CHECK(std::is_sorted(L.begin(), L.end()));
  // lower_bound also compares with a < b.
  auto iter = std::lower_bound(L.begin(), L.end(), get_vector<T>(1, 9));
  CHECK(iter != L.end() && *iter == get_vector<T>(1, 9));
}

// A pair of a matrix and a scalar as key: std::pair's operator< defers to the
// members, so the ordering of MyMatrix has to be reachable there too.
template <typename T> void test_set_pair() {
  std::set<std::pair<MyMatrix<T>, int>> S;
  S.insert({get_matrix<T>(1, 0, 0, 1), 1});
  S.insert({get_matrix<T>(1, 0, 0, 1), 1});
  S.insert({get_matrix<T>(1, 0, 0, 1), 2});
  CHECK(S.size() == 2);
}

template <typename T> void test_all(std::string const &name) {
  int before = n_fail;
  test_set_vector<T>();
  test_map_vector<T>();
  test_set_matrix<T>();
  test_sort_vector<T>();
  test_set_pair<T>();
  if (n_fail == before) {
    std::cerr << "Test_wasm_ordered_containers: " << name << " ok\n";
  }
}

int main() {
  test_all<cpp_rational>("cpp_rational");
  test_all<cpp_int>("cpp_int");
  if (n_fail == 0) {
    std::cerr << "Test_wasm_ordered_containers: all checks passed.\n";
    return 0;
  }
  std::cerr << "Test_wasm_ordered_containers: " << n_fail
            << " check(s) failed.\n";
  return 1;
}
