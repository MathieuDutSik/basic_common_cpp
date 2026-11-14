// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixModPrime.h"
#include "MAT_MatrixModInt64.h"
#include "Basic_functions.h"
#include <random>
#include <string>
#include <vector>
// clang-format on

/*
  Tests of MAT_MatrixModPrime.h and of the powers against independent
  computations: the rank of a square matrix against DeterminantMatModInt64,
  the reduced echelon shape, every input row recovered from its
  coordinates, the left nullspace by x M = 0 and its dimension, the
  solutions of x B = v for v in and out of the row space; ResPow, MyPow
  and MyPowSigned against repeated products, ResMul on any signs.
 */

int main() {
  try {
    std::mt19937_64 rng(1234);
    auto rnd = [&](int64_t lo, int64_t hi) -> int64_t {
      return lo + static_cast<int64_t>(rng() % static_cast<uint64_t>(hi - lo + 1));
    };
    int n_error = 0;
    auto fail = [&](std::string const &what) -> void {
      std::cerr << "ERROR: " << what << "\n";
      n_error++;
    };
    std::vector<int64_t> l_prime{2, 3, 5, 7, 101, 2147483629};
    for (int64_t p : l_prime) {
      for (int iter = 0; iter < 200; iter++) {
        int nr = static_cast<int>(rnd(1, 7));
        int nc = static_cast<int>(rnd(1, 7));
        // small entries make the rank deficiency frequent
        int64_t range = (iter % 2 == 0) ? 1 : 1000000;
        ModMatrix M(nr, std::vector<int64_t>(nc));
        for (auto &row : M) {
          for (auto &x : row) {
            x = rnd(-range, range);
          }
        }
        RowSpaceModP rs = ComputeRowSpaceModP(M, p);
        size_t rank = rs.pivots.size();
        // the echelon shape
        for (size_t i = 0; i < rank; i++) {
          for (size_t j = 0; j < rank; j++) {
            if (rs.basis[j][rs.pivots[i]] != ((i == j) ? 1 : 0)) {
              fail("echelon shape");
            }
          }
        }
        // every row from its coordinates
        for (auto const &row : M) {
          std::vector<int64_t> red(nc);
          for (int j = 0; j < nc; j++) {
            red[j] = ResMod(row[j], p);
          }
          std::vector<int64_t> c = CoordinatesInRowSpaceModP(rs, red);
          for (int j = 0; j < nc; j++) {
            int64_t s = 0;
            for (size_t i = 0; i < rank; i++) {
              s = ResAdd(s, ResMul(c[i], rs.basis[i][j], p), p);
            }
            if (s != red[j]) {
              fail("a row is not recovered from its coordinates");
            }
          }
        }
        // the rank of a square matrix against the determinant
        if (nr == nc && p < (int64_t(1) << 31)) {
          std::vector<int64_t> a;
          for (auto const &row : M) {
            for (auto x : row) {
              a.push_back(ResMod(x, p));
            }
          }
          bool full = DeterminantMatModInt64(a, nr, p, true) != 0;
          if (full != (static_cast<int>(rank) == nr)) {
            fail("the rank disagrees with the determinant");
          }
        }
        if (RankModPrime(M, p) != rank) {
          fail("RankModPrime");
        }
        // the left nullspace
        ModMatrix K = LeftNullspaceModP(M, p);
        if (K.size() != nr - rank || RankModPrime(K, p) != K.size()) {
          fail("the dimension of the left nullspace");
        }
        for (auto const &x : K) {
          for (int j = 0; j < nc; j++) {
            int64_t s = 0;
            for (int i = 0; i < nr; i++) {
              s = ResAdd(s, ResMul(x[i], ResMod(M[i][j], p), p), p);
            }
            if (s != 0) {
              fail("x M is not zero");
            }
          }
        }
        // x B = v on the echelon basis, v in the row space and maybe not
        if (rank > 0) {
          std::vector<int64_t> x0(rank), v(nc, 0);
          for (auto &x : x0) {
            x = rnd(0, p - 1);
          }
          for (size_t i = 0; i < rank; i++) {
            for (int j = 0; j < nc; j++) {
              v[j] = ResAdd(v[j], ResMul(x0[i], rs.basis[i][j], p), p);
            }
          }
          std::optional<std::vector<int64_t>> opt = SolveRowsModP(rs.basis, v, p);
          if (!opt || *opt != x0) {
            fail("SolveRowsModP on the row space");
          }
          if (static_cast<int>(rank) < nc) {
            std::vector<int64_t> w = v;
            for (int j = 0; j < nc; j++) {
              bool is_pivot = false;
              for (auto piv : rs.pivots) {
                if (static_cast<int>(piv) == j) {
                  is_pivot = true;
                }
              }
              if (!is_pivot) {
                w[j] = ResAdd(w[j], 1, p);
                break;
              }
            }
            if (SolveRowsModP(rs.basis, w, p)) {
              fail("SolveRowsModP off the row space");
            }
          }
        }
      }
    }
    // the powers
    for (int iter = 0; iter < 2000; iter++) {
      int64_t m = rnd(2, int64_t(1) << 40);
      int64_t a = rnd(-(int64_t(1) << 50), int64_t(1) << 50);
      int64_t e = rnd(0, 60);
      int64_t naive = ResMod(1, m);
      for (int64_t i = 0; i < e; i++) {
        naive = ResMul(naive, a, m);
      }
      if (ResPow(a, e, m) != naive) {
        fail("ResPow");
      }
      int64_t b = rnd(-(int64_t(1) << 50), int64_t(1) << 50);
      __int128 prod = static_cast<__int128>(a) * b;
      int64_t r = static_cast<int64_t>(((prod % m) + m) % m);
      if (ResMul(a, b, m) != r) {
        fail("ResMul with signs");
      }
    }
    for (int64_t a = -7; a <= 7; a++) {
      for (int e = 0; e <= 20; e++) {
        int64_t naive = 1;
        mpz_class naive_z(1);
        for (int i = 0; i < e; i++) {
          naive *= a;
          naive_z *= static_cast<long>(a);
        }
        if (MyPow<int64_t>(a, e) != naive || MyPow<mpz_class>(mpz_class(static_cast<long>(a)), e) != naive_z) {
          fail("MyPow");
        }
        if (a != 0) {
          mpq_class inv = mpq_class(1) / mpq_class(naive_z);
          if (MyPowSigned<mpq_class>(mpq_class(static_cast<long>(a)), -e) != inv) {
            fail("MyPowSigned");
          }
        }
      }
    }
    // a power that fits in int64 while the square after the last bit would not
    if (MyPow<int64_t>(3, 39) != 4052555153018976267LL) {
      fail("MyPow near the int64 bound");
    }
    std::cerr << "n_error=" << n_error << "\n";
    if (n_error > 0) {
      throw TerminalException{1};
    }
    std::cerr << "Normal termination of Test_MatrixModPrime\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_MatrixModPrime\n";
    exit(e.eVal);
  }
}
