// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// Correctness of the Eigen matrix products, in particular of those of
// fmpz_class / fmpq_class which go to flint (EigenFlintProduct.h): every
// form of product Eigen evaluates through its GEMM and GEMV kernels
// (transposed operands, blocks, row-major operands, strided vectors,
// destinations with a stride, +=, -=, scalar factors, chained products) is
// compared with a plain triple loop. Small and large sizes are covered,
// since below rows + cols + depth = 20 Eigen takes another path, and entries
// of ~90 bits, which flint multiplies by another algorithm. The other types
// check that the generic kernel is unchanged.
//
// Exits with a non-zero status on the first type with a failure.
// clang-format off
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "MAT_Matrix.h"
// clang-format on
#include <string>

namespace {

unsigned int test_seed = 17;

int next_entry(int spread) {
  test_seed = test_seed * 1103515245u + 12345u;
  return static_cast<int>((test_seed >> 16) % (2 * spread + 1)) - spread;
}

// Entries in [-7, 7], and for the field types denominators up to den.
template <typename T> T random_scalar(int den) {
  T x(next_entry(7));
  if constexpr (!is_implementation_of_Z<T>::value) {
    if (den > 1)
      x /= T(1 + (next_entry(den) + den) % den);
  }
  return x;
}

template <typename T> MyMatrix<T> random_matrix(int r, int c, int den) {
  MyMatrix<T> M(r, c);
  for (int i = 0; i < r; i++)
    for (int j = 0; j < c; j++)
      M(i, j) = random_scalar<T>(den);
  return M;
}

template <typename T> MyVector<T> random_vector(int n, int den) {
  MyVector<T> V(n);
  for (int i = 0; i < n; i++)
    V(i) = random_scalar<T>(den);
  return V;
}

// The reference: a triple loop on the coefficients, independent of the
// product kernels.
template <typename T>
MyMatrix<T> ref_product(MyMatrix<T> const &A, MyMatrix<T> const &B) {
  MyMatrix<T> C(A.rows(), B.cols());
  for (int i = 0; i < A.rows(); i++)
    for (int j = 0; j < B.cols(); j++) {
      T acc(0);
      for (int k = 0; k < A.cols(); k++)
        acc += A(i, k) * B(k, j);
      C(i, j) = acc;
    }
  return C;
}

template <typename T> MyMatrix<T> as_column(MyVector<T> const &V) {
  MyMatrix<T> M(V.size(), 1);
  for (int i = 0; i < V.size(); i++)
    M(i, 0) = V(i);
  return M;
}

struct Failures {
  std::string name;
  int n_fail = 0;
  int n_check = 0;
  template <typename T>
  void check(std::string const &what, int n, MyMatrix<T> const &X,
             MyMatrix<T> const &Y) {
    n_check++;
    if (X.rows() != Y.rows() || X.cols() != Y.cols() || X != Y) {
      std::cerr << "FAIL " << name << " " << what << " n=" << n << "\n";
      n_fail++;
    }
  }
};

template <typename T> void test_matrix_matrix(Failures &f, int n, int den) {
  using RowMat = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  MyMatrix<T> A = random_matrix<T>(n, n + 3, den);
  MyMatrix<T> B = random_matrix<T>(n + 3, n + 1, den);
  MyMatrix<T> C = random_matrix<T>(n + 1, n + 2, den);
  MyMatrix<T> AB = ref_product(A, B);
  MyMatrix<T> At = A.transpose(), Bt = B.transpose();
  f.check("A*B", n, MyMatrix<T>(A * B), AB);
  f.check("A*Bt^T", n, MyMatrix<T>(A * Bt.transpose()), AB);
  f.check("At^T*B", n, MyMatrix<T>(At.transpose() * B), AB);
  f.check("At^T*Bt^T", n, MyMatrix<T>(At.transpose() * Bt.transpose()), AB);
  f.check("(Bt*At)^T", n, MyMatrix<T>((Bt * At).transpose()), AB);
  f.check("A*B*C", n, MyMatrix<T>(A * B * C), ref_product(AB, C));
  f.check("A*(B*C)", n, MyMatrix<T>(A * (B * C)), ref_product(AB, C));
  RowMat Ar = A, Br = B;
  f.check("rowmajor*colmajor", n, MyMatrix<T>(Ar * B), AB);
  f.check("colmajor*rowmajor", n, MyMatrix<T>(A * Br), AB);
  RowMat ABr = Ar * Br;
  f.check("rowmajor result", n, MyMatrix<T>(ABr), AB);
  // Accumulation and scalar factors.
  MyMatrix<T> D = AB;
  D += A * B;
  f.check("+=", n, D, MyMatrix<T>(T(2) * AB));
  D -= A * B;
  f.check("-=", n, D, AB);
  D.noalias() -= T(3) * (A * B);
  f.check("noalias -= 3*", n, D, MyMatrix<T>(T(-2) * AB));
  f.check("3*(A*B)", n, MyMatrix<T>(T(3) * (A * B)), MyMatrix<T>(T(3) * AB));
  f.check("(3*A)*B", n, MyMatrix<T>((T(3) * A) * B), MyMatrix<T>(T(3) * AB));
  // Blocks as operands and as destination (a stride larger than the rows).
  if (n >= 2) {
    MyMatrix<T> Ab = A.block(1, 1, n - 1, n), Bb = B.block(1, 0, n, n);
    f.check("block*block", n, MyMatrix<T>(A.block(1, 1, n - 1, n) * B.block(1, 0, n, n)),
            ref_product(Ab, Bb));
    f.check("topRows*B", n, MyMatrix<T>(A.topRows(n - 1) * B),
            ref_product(MyMatrix<T>(A.topRows(n - 1)), B));
    MyMatrix<T> E = MyMatrix<T>::Zero(n + 2, n + 3);
    E.block(1, 2, n, n + 1).noalias() = A * B;
    f.check("block destination", n, MyMatrix<T>(E.block(1, 2, n, n + 1)), AB);
    f.check("block destination untouched", n,
            MyMatrix<T>(E.topRows(1)), MyMatrix<T>(MyMatrix<T>::Zero(1, n + 3)));
    E.block(1, 2, n, n + 1).noalias() += A * B;
    f.check("block destination +=", n, MyMatrix<T>(E.block(1, 2, n, n + 1)),
            MyMatrix<T>(T(2) * AB));
  }
  // Aliasing: the product into one of its operands.
  MyMatrix<T> S = random_matrix<T>(n, n, den);
  MyMatrix<T> S2 = ref_product(S, S);
  S = S * S;
  f.check("S = S*S", n, S, S2);
}

template <typename T> void test_matrix_vector(Failures &f, int n, int den) {
  MyMatrix<T> A = random_matrix<T>(n, n + 3, den);
  MyVector<T> v = random_vector<T>(n + 3, den);
  MyVector<T> w = random_vector<T>(n, den);
  MyMatrix<T> Av = ref_product(A, as_column(v));
  MyMatrix<T> Atw = ref_product(MyMatrix<T>(A.transpose()), as_column(w));
  f.check("A*v", n, as_column<T>(A * v), Av);
  f.check("w^T*A", n, MyMatrix<T>(MyVector<T>((w.transpose() * A).transpose())), Atw);
  f.check("A^T*w", n, as_column<T>(A.transpose() * w), Atw);
  MyMatrix<T> M = random_matrix<T>(n + 3, 4, den);
  f.check("A*col", n, as_column<T>(A * M.col(2)),
          ref_product(A, MyMatrix<T>(M.col(2))));
  // A row of a column-major matrix: entries with a stride.
  MyMatrix<T> Mt = M.transpose();
  f.check("A*row^T (strided)", n, as_column<T>(A * Mt.row(2).transpose()),
          ref_product(A, MyMatrix<T>(M.col(2))));
  MyVector<T> u = A * v;
  u += A * v;
  f.check("v +=", n, as_column(u), MyMatrix<T>(T(2) * Av));
  u -= T(5) * (A * v);
  f.check("v -= 5*", n, as_column(u), MyMatrix<T>(T(-3) * Av));
  // The result into a row of a matrix: a destination with a stride.
  MyMatrix<T> D = MyMatrix<T>::Zero(3, n + 3);
  D.row(1).noalias() = w.transpose() * A;
  f.check("row destination", n, MyMatrix<T>(D.row(1).transpose()), Atw);
  D.row(1).noalias() += w.transpose() * A;
  f.check("row destination +=", n, MyMatrix<T>(D.row(1).transpose()),
          MyMatrix<T>(T(2) * Atw));
  f.check("row destination untouched", n, MyMatrix<T>(D.row(0)),
          MyMatrix<T>(MyMatrix<T>::Zero(1, n + 3)));
}

template <typename T> void test_large_entries(Failures &f, int den) {
  // Entries of about 90 bits, a few of them negative.
  for (int n : {3, 12, 25}) {
    MyMatrix<T> A = random_matrix<T>(n, n, den);
    MyMatrix<T> B = random_matrix<T>(n, n, den);
    T big = T(1000000007) * T(1000000009) * T(998244353);
    for (int i = 0; i < n; i++) {
      A(i, (i + 1) % n) *= big;
      B(i, i) = B(i, i) * big + T(1);
    }
    f.check("big A*B", n, MyMatrix<T>(A * B), ref_product(A, B));
    MyVector<T> v = GetMatrixRow(B, 0);
    f.check("big A*v", n, as_column<T>(A * v), ref_product(A, as_column(v)));
  }
}

template <typename T> int test_type(std::string const &name, int den) {
  Failures f;
  f.name = name;
  for (int n : {0, 1, 2, 3, 5, 8, 13, 19, 31, 50}) {
    test_matrix_matrix<T>(f, n, den);
    test_matrix_vector<T>(f, n, den);
  }
  if constexpr (!std::is_integral_v<T>)
    test_large_entries<T>(f, den);
  std::cerr << name << ": " << f.n_check << " checks, " << f.n_fail
            << " failures\n";
  return f.n_fail;
}

}  // namespace

int main() {
  HumanTime time;
  try {
    int n_fail = 0;
#ifdef ENABLE_FLINT_SUPPORT
    n_fail += test_type<fmpz_class>("fmpz_class", 1);
    n_fail += test_type<fmpq_class>("fmpq_class", 1);
    n_fail += test_type<fmpq_class>("fmpq_class(den)", 6);
#else
    std::cerr << "Built without ENABLE_FLINT_SUPPORT: the flint types are not tested\n";
#endif
    n_fail += test_type<mpz_class>("mpz_class", 1);
    n_fail += test_type<mpq_class>("mpq_class(den)", 6);
    n_fail += test_type<int64_t>("int64_t", 1);
    if (n_fail > 0) {
      std::cerr << "Test_EigenFlintProduct: " << n_fail << " failures\n";
      return 1;
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
