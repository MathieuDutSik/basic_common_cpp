// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// Benchmark of the matrix kernels across the exact arithmetics, in
// particular gmp (mpz_class / mpq_class) against flint (fmpz_class /
// fmpq_class, with ENABLE_FLINT_SUPPORT). One arithmetic per invocation so
// that the allocator state of one type cannot influence the timings of
// another; the inputs are rebuilt identically in every invocation from a
// fixed seed, so every arithmetic computes on the same matrices.
//
// The operations, chosen to exercise the two regimes of an exact
// computation:
//   product       entries stay small          (the small-integer regime)
//   product_big   ~90 bit entries             (the large-integer regime)
//   determinant   Bareiss, coefficient growth
//   inverse       field types only
//   nullspace     field types only
//   hnf           integer rings only
//
// Each operation is run several times and the minimum is reported, which is
// the standard way of removing scheduling noise from a benchmark.
//
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "MAT_MatrixInt.h"
// clang-format on
#include <chrono>
#include <random>
#include <string>
#include <vector>

static auto now() { return std::chrono::steady_clock::now(); }
static double ms(std::chrono::steady_clock::duration d) {
  return std::chrono::duration<double, std::milli>(d).count();
}

std::vector<int> random_entries(size_t cnt, int spread, unsigned int seed) {
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int> dist(-spread, spread);
  std::vector<int> ent(cnt);
  for (auto &x : ent)
    x = dist(gen);
  return ent;
}

template <typename T>
MyMatrix<T> build_matrix(int rows, int cols, std::vector<int> const &ent) {
  MyMatrix<T> M(rows, cols);
  size_t pos = 0;
  for (int i = 0; i < rows; i++)
    for (int j = 0; j < cols; j++)
      M(i, j) = ent[pos++];
  return M;
}

// A rank r matrix as a product of two random small matrices, the standard
// way of getting a non-trivial nullspace.
template <typename T>
MyMatrix<T> build_rank_deficient(int rows, int cols, int r,
                                 unsigned int seed) {
  MyMatrix<T> B =
      build_matrix<T>(rows, r, random_entries(rows * r, 5, seed));
  MyMatrix<T> C =
      build_matrix<T>(r, cols, random_entries(r * cols, 5, seed + 1));
  return B * C;
}

template <typename F>
double time_best(int reps, F f) {
  double best = -1;
  for (int i = 0; i < reps; i++) {
    auto t0 = now();
    f();
    double dur = ms(now() - t0);
    if (best < 0 || dur < best)
      best = dur;
  }
  return best;
}

template <typename T> void bench_arith(std::string const &name) {
  auto report = [&](std::string const &op, double dur) {
    std::cout << name << " " << op << " " << dur << " ms" << std::endl;
  };
  auto announce = [&](std::string const &op) {
    std::cerr << "start " << op << "\n";
  };
  // product, small entries
  {
    int n = 120;
    MyMatrix<T> A = build_matrix<T>(n, n, random_entries(n * n, 100, 71));
    MyMatrix<T> B = build_matrix<T>(n, n, random_entries(n * n, 100, 72));
    MyMatrix<T> C;
    announce("product");
    report("product", time_best(5, [&]() { C = A * B; }));
  }
  // product, ~90 bit entries
  {
    int n = 80;
    std::vector<int> e1 = random_entries(n * n, 1000000000, 73);
    std::vector<int> e2 = random_entries(n * n, 1000000000, 74);
    std::vector<int> e3 = random_entries(n * n, 1000000000, 75);
    auto build_big = [&](std::vector<int> const &u, std::vector<int> const &v,
                         std::vector<int> const &w) {
      MyMatrix<T> M(n, n);
      size_t pos = 0;
      for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
          M(i, j) = T(u[pos]) * T(v[pos]) * T(w[pos]);
          pos++;
        }
      return M;
    };
    MyMatrix<T> A = build_big(e1, e2, e3);
    MyMatrix<T> B = build_big(e2, e3, e1);
    MyMatrix<T> C;
    announce("product_big");
    report("product_big", time_best(5, [&]() { C = A * B; }));
  }
  // determinant
  {
    int n = 60;
    MyMatrix<T> A = build_matrix<T>(n, n, random_entries(n * n, 10, 76));
    T det;
    announce("determinant");
    report("determinant", time_best(3, [&]() { det = DeterminantMat(A); }));
  }
  if constexpr (is_ring_field<T>::value) {
    // inverse
    {
      int n = 50;
      MyMatrix<T> A = build_matrix<T>(n, n, random_entries(n * n, 10, 77));
      MyMatrix<T> Inv;
      announce("inverse");
      report("inverse", time_best(3, [&]() { Inv = Inverse(A); }));
    }
    // nullspace
    {
      MyMatrix<T> A = build_rank_deficient<T>(60, 90, 40, 78);
      MyMatrix<T> Ker;
      announce("nullspace");
      report("nullspace", time_best(3, [&]() { Ker = NullspaceMat(A); }));
    }
  } else {
    // hnf: the coefficient growth of the transformation matrix explodes
    // past n ~ 35, so the size stays below that.
    {
      int n = 32;
      MyMatrix<T> A = build_matrix<T>(n, n, random_entries(n * n, 10, 79));
      announce("hnf");
      report("hnf", time_best(3, [&]() {
               std::pair<MyMatrix<T>, MyMatrix<T>> ePair =
                   ComputeRowHermiteNormalForm(A);
             }));
    }
  }
}

void process(std::string const &arith) {
#ifndef DISABLE_GMP_ARITHMETIC
  if (arith == "integer") {
    return bench_arith<mpz_class>(arith);
  }
  if (arith == "rational") {
    return bench_arith<mpq_class>(arith);
  }
#endif
#ifdef ENABLE_FLINT_SUPPORT
  if (arith == "flint_integer") {
    return bench_arith<fmpz_class>(arith);
  }
  if (arith == "flint_rational") {
    return bench_arith<fmpq_class>(arith);
  }
#endif
  std::cerr << "Failed to find a matching entry for arith\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  try {
    if (argc != 2) {
      std::cerr << "This program is used as\n";
      std::cerr << "Bench_matrix_arithmetic [arith]\n";
      std::cerr << "    where\n";
      std::cerr << "arith: integer, rational";
      std::cerr << ", flint_integer, flint_rational (with flint support)\n";
      return -1;
    }
    std::string arith = argv[1];
    process(arith);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
}
