// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheory.h"
#include "NumberTheorySafeInt.h"
#include "MAT_MatrixMod.h"
#include <random>
#include <string>
#include <vector>
// clang-format on

/*
  Tests of RecSolutionMatMod, the membership in a space given modulo p.

  Every test below decides the same question twice: once with the
  structure, and once with an independent reference. The reference is the
  rank, computed by SelectRowColMatMod: a vector V belongs to the space
  spanned by the rows of Space if and only if adding V as a row does not
  raise the rank. That reference shares the elimination with the
  structure but nothing of the reformulation by equations, which is what
  the tests are about.
 */

static std::mt19937 &get_generator() {
  // Fixed seed: a failure has to be reproducible, and the point of the
  // test is coverage over many shapes rather than a new draw each run.
  static std::mt19937 gen(20250919);
  return gen;
}

template <typename T>
MyMatrix<T> get_random_matrix(int n_row, int n_col, T const &TheMod) {
  int mod_i = UniversalScalarConversion<int, T>(TheMod);
  std::uniform_int_distribution<int> distr(0, mod_i - 1);
  MyMatrix<T> M(n_row, n_col);
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n_col; j++) {
      M(i, j) = T(distr(get_generator()));
    }
  }
  return M;
}

template <typename T> MyVector<T> get_random_vector(int n, T const &TheMod) {
  MyMatrix<T> M = get_random_matrix<T>(1, n, TheMod);
  return GetMatrixRow(M, 0);
}

// The reference answer: V is in the span of the rows of Space if and
// only if it does not increase the rank.
template <typename T>
bool reference_membership(MyMatrix<T> const &Space, MyVector<T> const &V,
                          T const &TheMod) {
  int n_row = Space.rows();
  int n = Space.cols();
  MyMatrix<T> Ext(n_row + 1, n);
  for (int i = 0; i < n_row; i++) {
    for (int j = 0; j < n; j++) {
      Ext(i, j) = Space(i, j);
    }
  }
  for (int j = 0; j < n; j++) {
    Ext(n_row, j) = V(j);
  }
  size_t rnk_space = SelectRowColMatMod(Space, TheMod).TheRank;
  size_t rnk_ext = SelectRowColMatMod(Ext, TheMod).TheRank;
  return rnk_space == rnk_ext;
}

static void report_failure(std::string const &test, std::string const &detail) {
  std::cerr << "TESTMOD: FAILURE in " << test << ": " << detail << "\n";
  throw TerminalException{1};
}

/*
  A combination of the rows of the space is in the space. This is the
  direction that a structure returning always false would pass, so it is
  checked together with the reference below rather than on its own.
 */
template <typename T>
void test_combinations_are_contained(MyMatrix<T> const &Space, T const &TheMod,
                                     int n_test) {
  int n_row = Space.rows();
  int n = Space.cols();
  RecSolutionMatMod<T> rec(Space, TheMod);
  for (int i_test = 0; i_test < n_test; i_test++) {
    MyVector<T> coeff = get_random_vector<T>(n_row, TheMod);
    MyVector<T> V = ZeroVector<T>(n);
    for (int i = 0; i < n_row; i++) {
      for (int j = 0; j < n; j++) {
        T val = V(j) + coeff(i) * Space(i, j);
        V(j) = ResInt(val, TheMod);
      }
    }
    if (!rec.has_solution_v(V)) {
      report_failure("test_combinations_are_contained",
                     "a combination of the rows was rejected");
    }
  }
}

/*
  The structure and the rank reference agree on vectors drawn at random.
  With a space of dimension d inside (Z/pZ)^n, such a vector misses the
  space with probability 1 - p^(d-n), so both answers occur over a run;
  the two are counted so that the test cannot quietly degenerate into
  checking one of the two directions only.
 */
static size_t n_random_inside = 0;
static size_t n_random_outside = 0;

template <typename T>
void test_against_reference(MyMatrix<T> const &Space, T const &TheMod,
                            int n_test) {
  int n = Space.cols();
  RecSolutionMatMod<T> rec(Space, TheMod);
  for (int i_test = 0; i_test < n_test; i_test++) {
    MyVector<T> V = get_random_vector<T>(n, TheMod);
    bool obtained = rec.has_solution_v(V);
    bool expected = reference_membership(Space, V, TheMod);
    if (obtained != expected) {
      report_failure("test_against_reference",
                     "the structure and the rank disagree on a random vector");
    }
    if (expected) {
      n_random_inside++;
    } else {
      n_random_outside++;
    }
  }
}

/*
  The space contains itself, and is_containing_m is the conjunction of
  has_solution_v over the rows.
 */
template <typename T>
void test_is_containing_m(MyMatrix<T> const &Space, T const &TheMod,
                          int n_test) {
  int n_row = Space.rows();
  int n = Space.cols();
  RecSolutionMatMod<T> rec(Space, TheMod);
  if (n_row > 0 && !rec.is_containing_m(Space)) {
    report_failure("test_is_containing_m", "the space does not contain itself");
  }
  for (int i_test = 0; i_test < n_test; i_test++) {
    int n_row_test = 1 + (i_test % 3);
    MyMatrix<T> M = get_random_matrix<T>(n_row_test, n, TheMod);
    bool obtained = rec.is_containing_m(M);
    bool expected = true;
    for (int i = 0; i < n_row_test; i++) {
      MyVector<T> V = GetMatrixRow(M, i);
      if (!rec.has_solution_v(V)) {
        expected = false;
      }
    }
    if (obtained != expected) {
      report_failure("test_is_containing_m",
                     "is_containing_m is not the conjunction over the rows");
    }
  }
}

/*
  The answer depends on the vector modulo p only, so shifting the entries
  by multiples of p may not change it.

  The multiples are taken large on purpose. Whether the structure reduces
  its argument or not, the answer comes out the same, since the terms it
  adds are multiples of p and are reduced away at the end; what the
  reduction buys is the size of what is handled on the way. Without it a
  shifted entry of about 10^17 meets a coefficient of up to p - 1 and the
  product leaves the range of a 64 bit integer, so on SafeInt64 a missing
  reduction surfaces here as an overflow. On the unbounded types this
  test only checks the answer.
 */
template <typename T>
void test_reduction_is_applied(MyMatrix<T> const &Space, T const &TheMod,
                               int n_test) {
  int n = Space.cols();
  RecSolutionMatMod<T> rec(Space, TheMod);
  int64_t max_shift = 1000000000000000;
  std::uniform_int_distribution<int64_t> distr(-max_shift, max_shift);
  for (int i_test = 0; i_test < n_test; i_test++) {
    MyVector<T> V = get_random_vector<T>(n, TheMod);
    MyVector<T> Vshift(n);
    for (int i = 0; i < n; i++) {
      int64_t shift = distr(get_generator());
      T shift_T = UniversalScalarConversion<T, int64_t>(shift);
      Vshift(i) = V(i) + shift_T * TheMod;
    }
    if (rec.has_solution_v(V) != rec.has_solution_v(Vshift)) {
      report_failure("test_reduction_is_applied",
                     "shifting a vector by multiples of the modulus changed "
                     "the answer");
    }
  }
}

/*
  The two extreme spaces. A space of n independent rows is the whole of
  (Z/pZ)^n and contains everything; a space of no rows at all is the zero
  space and contains only the zero vector.
 */
template <typename T> void test_extreme_spaces(int n, T const &TheMod) {
  MyMatrix<T> Full = IdentityMat<T>(n);
  RecSolutionMatMod<T> rec_full(Full, TheMod);
  for (int i_test = 0; i_test < 10; i_test++) {
    MyVector<T> V = get_random_vector<T>(n, TheMod);
    if (!rec_full.has_solution_v(V)) {
      report_failure("test_extreme_spaces", "the full space rejected a vector");
    }
  }
  MyMatrix<T> Empty(0, n);
  RecSolutionMatMod<T> rec_empty(Empty, TheMod);
  MyVector<T> Zero = ZeroVector<T>(n);
  if (!rec_empty.has_solution_v(Zero)) {
    report_failure("test_extreme_spaces",
                   "the zero space rejected the zero vector");
  }
  for (int i = 0; i < n; i++) {
    MyVector<T> V = ZeroVector<T>(n);
    V(i) = T(1);
    if (rec_empty.has_solution_v(V)) {
      report_failure("test_extreme_spaces",
                     "the zero space accepted a nonzero vector");
    }
  }
}

/*
  A space built with a prescribed number of independent rows, so that the
  cases where the spanning family is not free are covered as well: the
  rows beyond the first n_gen are combinations of them.
 */
template <typename T>
MyMatrix<T> get_space_with_redundancy(int n_row, int n_gen, int n,
                                      T const &TheMod) {
  MyMatrix<T> Space = get_random_matrix<T>(n_row, n, TheMod);
  for (int i = n_gen; i < n_row; i++) {
    for (int j = 0; j < n; j++) {
      Space(i, j) = 0;
    }
    for (int i_gen = 0; i_gen < n_gen; i_gen++) {
      MyVector<T> coeff = get_random_vector<T>(1, TheMod);
      for (int j = 0; j < n; j++) {
        T val = Space(i, j) + coeff(0) * Space(i_gen, j);
        Space(i, j) = ResInt(val, TheMod);
      }
    }
  }
  return Space;
}

template <typename T> void full_process_type(int n_iter) {
  std::vector<int> ListMod{2, 3, 5, 7, 11, 101};
  int n_test = 20;
  for (auto &mod_i : ListMod) {
    T TheMod(mod_i);
    std::cerr << "TESTMOD: TheMod=" << TheMod << "\n";
    test_extreme_spaces<T>(4, TheMod);
    for (int i_iter = 0; i_iter < n_iter; i_iter++) {
      // The dimensions cycle over the shapes that matter: fewer rows
      // than columns, as many, and more.
      int n = 1 + (i_iter % 6);
      int n_row = 1 + ((i_iter / 6) % 7);
      int n_gen = 1 + (i_iter % n_row);
      MyMatrix<T> Space = get_space_with_redundancy<T>(n_row, n_gen, n, TheMod);
      test_combinations_are_contained(Space, TheMod, n_test);
      test_against_reference(Space, TheMod, n_test);
      test_is_containing_m(Space, TheMod, n_test);
      test_reduction_is_applied(Space, TheMod, n_test);
    }
  }
  std::cerr << "TESTMOD: random vectors inside=" << n_random_inside
            << " outside=" << n_random_outside << "\n";
  if (n_random_inside == 0 || n_random_outside == 0) {
    report_failure("full_process_type",
                   "the random vectors all fell on the same side, so the "
                   "comparison with the reference tested one direction only");
  }
  std::cerr << "TESTMOD: all the RecSolutionMatMod tests passed\n";
}

void process(std::string const &arith, int n_iter) {
  if (arith == "safe_integer") {
    using T = SafeInt64;
    return full_process_type<T>(n_iter);
  }
  if (arith == "mpz_class") {
    using T = mpz_class;
    return full_process_type<T>(n_iter);
  }
  if (arith == "boost_cpp_int") {
    using T = boost::multiprecision::cpp_int;
    return full_process_type<T>(n_iter);
  }
  std::cerr << "Failed to find a matching entry for arith\n";
  std::cerr << "Allowed values are safe_integer, mpz_class, boost_cpp_int\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 2 && argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "Test_MatrixMod [arith]\n";
      std::cerr << "or\n";
      std::cerr << "Test_MatrixMod [arith] [n_iter]\n";
      std::cerr << "---\n";
      std::cerr << "arith  : safe_integer, mpz_class, boost_cpp_int\n";
      std::cerr << "n_iter : number of spaces tested per modulus "
                << "(default 42)\n";
      return -1;
    }
    std::string arith = argv[1];
    int n_iter = 42;
    if (argc == 3) {
      n_iter = std::stoi(argv[2]);
    }
    process(arith, n_iter);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened in the computation\n";
    exit(e.eVal);
  }
  runtime(time);
}
