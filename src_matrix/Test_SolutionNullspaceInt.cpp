// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheory.h"
#include "NumberTheorySafeInt.h"
#include "MAT_MatrixInt.h"
// clang-format on

// Tests of the linear system solving and the integral kernels:
// --- SolutionMat / SolutionIntMat: a planted solution is recovered
//     (verified by substitution), an unsolvable system is rejected, and
//     over a ring SolutionMat returns the fractional solution over the
//     overlying field whenever one exists.
// --- NullspaceIntTrMat / NullspaceIntMat: the kernel vectors annihilate
//     the matrix, the rank is the nullity, and the returned basis is
//     saturated (every integral kernel vector is an integral combination
//     of the basis).
// --- GetZbasis: double inclusion between the generators and the basis.
// --- IntersectionLattice: the intersection is inside both lattices and
//     contains det(M2) L1 and det(M1) L2.

template <typename T> MyMatrix<T> RandomMatrix(int n_row, int n_col, int amp) {
  MyMatrix<T> A(n_row, n_col);
  for (int i = 0; i < n_row; i++)
    for (int j = 0; j < n_col; j++)
      A(i, j) = T((random() % (2 * amp + 1)) - amp);
  return A;
}

template <typename T> MyVector<T> RandomVector(int n, int amp) {
  MyVector<T> V(n);
  for (int i = 0; i < n; i++)
    V(i) = T((random() % (2 * amp + 1)) - amp);
  return V;
}

template <typename T> void process_solution(int n, int nb) {
  for (int i = 0; i < nb; i++) {
    int n_row = n + (i % 3) - 1;
    MyMatrix<T> A = RandomMatrix<T>(n_row, n, 5);
    // A planted integral solution x A = b has to be recovered.
    MyVector<T> x = RandomVector<T>(n_row, 4);
    MyVector<T> b = A.transpose() * x;
    std::optional<MyVector<T>> opt_int = SolutionIntMat(A, b);
    if (!opt_int) {
      std::cerr << "SolutionIntMat missed a solvable system\n";
      std::cerr << "A=\n";
      WriteMatrix(std::cerr, A);
      throw TerminalException{1};
    }
    MyVector<T> b_check = A.transpose() * (*opt_int);
    if (b_check != b) {
      std::cerr << "The solution of SolutionIntMat is incorrect\n";
      std::cerr << "A=\n";
      WriteMatrix(std::cerr, A);
      throw TerminalException{1};
    }
    // Over the ring, SolutionMat returns the (fractional in general)
    // solution over the overlying field, and misses nothing solvable.
    using Tfield = typename overlying_field<T>::field_type;
    MyMatrix<Tfield> Af = UniversalMatrixConversion<Tfield, T>(A);
    MyVector<Tfield> bf = UniversalVectorConversion<Tfield, T>(b);
    std::optional<MyVector<Tfield>> opt_ring = SolutionMat(A, b);
    if (!opt_ring) {
      std::cerr << "SolutionMat missed a solvable system\n";
      throw TerminalException{1};
    }
    MyVector<Tfield> b_check_fld = Af.transpose() * (*opt_ring);
    if (b_check_fld != bf) {
      std::cerr << "The solution of SolutionMat is incorrect\n";
      throw TerminalException{1};
    }
    // A right hand side outside of the row space has to be rejected.
    if (RankMat(A) < n) {
      MyVector<T> b_bad = RandomVector<T>(n, 5);
      MyMatrix<T> Aext(n_row + 1, n);
      for (int j = 0; j < n; j++) {
        for (int k = 0; k < n_row; k++)
          Aext(k, j) = A(k, j);
        Aext(n_row, j) = b_bad(j);
      }
      if (RankMat(Aext) > RankMat(A)) {
        if (SolutionMat(A, b_bad)) {
          std::cerr << "SolutionMat accepted an unsolvable system\n";
          throw TerminalException{1};
        }
        if (SolutionIntMat(A, b_bad)) {
          std::cerr << "SolutionIntMat accepted an unsolvable system\n";
          throw TerminalException{1};
        }
      }
    }
  }
  // A rationally solvable but integrally unsolvable system.
  MyMatrix<T> Adouble = 2 * IdentityMat<T>(n);
  MyVector<T> e1 = ZeroVector<T>(n);
  e1(0) = T(1);
  if (SolutionIntMat(Adouble, e1)) {
    std::cerr << "SolutionIntMat accepted a non-integral system\n";
    throw TerminalException{1};
  }
  // The ring SolutionMat has to return the fractional solution.
  using Tfield = typename overlying_field<T>::field_type;
  std::optional<MyVector<Tfield>> opt_half = SolutionMat(Adouble, e1);
  if (!opt_half) {
    std::cerr << "SolutionMat rejected a rationally solvable system\n";
    throw TerminalException{1};
  }
}

template <typename T> void process_nullspace_int(int n, int nb) {
  for (int i = 0; i < nb; i++) {
    // A matrix of controlled rank r in dimension n x n_col.
    int n_col = n + 1 + (i % 2);
    int r = 1 + (random() % (n - 1));
    MyMatrix<T> B = RandomMatrix<T>(r, n_col, 4);
    MyMatrix<T> C = RandomMatrix<T>(n, r, 3);
    MyMatrix<T> M = C * B;
    int rank = RankMat(M);
    // The rows v of NullspaceIntMat satisfy v M = 0.
    MyMatrix<T> NSProw = NullspaceIntMat(M);
    if (NSProw.rows() != M.rows() - rank) {
      std::cerr << "NullspaceIntMat has the wrong dimension\n";
      throw TerminalException{1};
    }
    for (int k = 0; k < NSProw.rows(); k++) {
      MyVector<T> v = GetMatrixRow(NSProw, k);
      MyVector<T> prod = M.transpose() * v;
      if (!IsZeroVector(prod)) {
        std::cerr << "A row of NullspaceIntMat is not in the kernel\n";
        throw TerminalException{1};
      }
    }
    // The rows v of NullspaceIntTrMat satisfy M v = 0.
    MyMatrix<T> NSPcol = NullspaceIntTrMat(M);
    if (NSPcol.rows() != M.cols() - rank) {
      std::cerr << "NullspaceIntTrMat has the wrong dimension\n";
      throw TerminalException{1};
    }
    for (int k = 0; k < NSPcol.rows(); k++) {
      MyVector<T> v = GetMatrixRow(NSPcol, k);
      MyVector<T> prod = M * v;
      if (!IsZeroVector(prod)) {
        std::cerr << "A row of NullspaceIntTrMat is not in the kernel\n";
        throw TerminalException{1};
      }
    }
    // Saturation: an integral kernel vector is an integral combination.
    if (NSProw.rows() > 0) {
      MyVector<T> w = ZeroVector<T>(M.rows());
      for (int k = 0; k < NSProw.rows(); k++) {
        T c = T((random() % 7) - 3);
        for (int j = 0; j < M.rows(); j++)
          w(j) += c * NSProw(k, j);
      }
      if (!SolutionIntMat(NSProw, w)) {
        std::cerr << "The saturation of NullspaceIntMat fails\n";
        throw TerminalException{1};
      }
    }
  }
}

template <typename T> void process_zbasis(int n, int nb) {
  for (int i = 0; i < nb; i++) {
    // Generators: r independent rows and integer combinations of them.
    int r = 1 + (random() % n);
    int m = r + 2 + (random() % 3);
    MyMatrix<T> B = RandomMatrix<T>(r, n, 4);
    MyMatrix<T> Gens(m, n);
    for (int k = 0; k < m; k++) {
      if (k < r) {
        Gens.row(k) = B.row(k);
      } else {
        MyVector<T> comb = ZeroVector<T>(n);
        for (int l = 0; l < r; l++) {
          T c = T((random() % 5) - 2);
          for (int j = 0; j < n; j++)
            comb(j) += c * B(l, j);
        }
        for (int j = 0; j < n; j++)
          Gens(k, j) = comb(j);
      }
    }
    MyMatrix<T> Basis = GetZbasis(Gens);
    if (RankMat(Basis) != RankMat(Gens) || Basis.rows() != RankMat(Gens)) {
      std::cerr << "GetZbasis has the wrong rank\n";
      throw TerminalException{1};
    }
    for (int k = 0; k < m; k++) {
      MyVector<T> v = GetMatrixRow(Gens, k);
      if (!SolutionIntMat(Basis, v)) {
        std::cerr << "A generator is not in the Z-span of the basis\n";
        throw TerminalException{1};
      }
    }
    for (int k = 0; k < Basis.rows(); k++) {
      MyVector<T> v = GetMatrixRow(Basis, k);
      if (!SolutionIntMat(Gens, v)) {
        std::cerr << "A basis vector is not in the Z-span of the generators\n";
        throw TerminalException{1};
      }
    }
  }
}

template <typename T> void process_intersection(int n, int nb) {
  for (int i = 0; i < nb; i++) {
    MyMatrix<T> M1(n, n), M2(n, n);
    while (true) {
      M1 = RandomMatrix<T>(n, n, 4);
      if (DeterminantMat(M1) != 0)
        break;
    }
    while (true) {
      M2 = RandomMatrix<T>(n, n, 4);
      if (DeterminantMat(M2) != 0)
        break;
    }
    MyMatrix<T> L = IntersectionLattice(M1, M2);
    if (RankMat(L) != n) {
      std::cerr << "The intersection lattice has the wrong rank\n";
      throw TerminalException{1};
    }
    // The intersection is inside both lattices.
    for (int k = 0; k < L.rows(); k++) {
      MyVector<T> v = GetMatrixRow(L, k);
      if (!SolutionIntMat(M1, v) || !SolutionIntMat(M2, v)) {
        std::cerr << "The intersection is not in both lattices\n";
        throw TerminalException{1};
      }
    }
    // det(M2) L1 is inside both lattices, hence in the intersection.
    T det2 = DeterminantMat(M2);
    MyVector<T> x = RandomVector<T>(n, 3);
    MyVector<T> v = det2 * (M1.transpose() * x);
    if (!SolutionIntMat(L, v)) {
      std::cerr << "The intersection lattice misses a common vector\n";
      throw TerminalException{1};
    }
  }
}

template <typename T> void process(int n) {
  int nb = 30;
  process_solution<T>(n, nb);
  std::cerr << "process_solution done\n";
  process_nullspace_int<T>(n, nb);
  std::cerr << "process_nullspace_int done\n";
  process_zbasis<T>(n, nb);
  std::cerr << "process_zbasis done\n";
  process_intersection<T>(n, nb);
  std::cerr << "process_intersection done\n";
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "Test_SolutionNullspaceInt [arith] [n]\n";
      return -1;
    }
    std::string arith = argv[1];
    int n;
    sscanf(argv[2], "%d", &n);
    auto f = [&]() -> void {
      if (arith == "mpz_class")
        return process<mpz_class>(n);
      if (arith == "safe_integer")
        return process<SafeInt64>(n);
      if (arith == "boost_cpp_int")
        return process<boost::multiprecision::cpp_int>(n);
      std::cerr << "Failed to find a matching entry\n";
      throw TerminalException{1};
    };
    f();
    std::cerr << "Normal termination of Test_SolutionNullspaceInt\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of Test_SolutionNullspaceInt\n";
    exit(e.eVal);
  }
  runtime(time);
}
