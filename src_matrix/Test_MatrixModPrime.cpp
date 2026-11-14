// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixModPrime.h"
#include "MAT_MatrixModInt64.h"
#include "Basic_functions.h"
#include <numeric>
#include <random>
#include <string>
#include <vector>
// clang-format on

/*
  Tests of MAT_MatrixModPrime.h, ResidueArithmetic.h and of the powers
  against independent computations: the rank of a square matrix against
  DeterminantMatModInt64, the reduced echelon shape, every input row
  recovered from its coordinates, the free kernel basis and the left
  nullspace by their products with M and their dimensions, the solutions
  of x B = v for v in and out of the row space; ResModSigned, ResModOf,
  ResMul, ResAdd, ResSub, ResPow, ResInv, MyPow and MyPowSigned against
  __int128, mpz_class and repeated products.
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
    std::vector<T_small_prime> l_prime{2, 3, 5, 7, 101, 2147483629,
                                       (T_small_prime(1) << 61) - 1};
    for (T_small_prime p : l_prime) {
      for (int iter = 0; iter < 200; iter++) {
        int nr = static_cast<int>(rnd(1, 7));
        int nc = static_cast<int>(rnd(1, 7));
        // small entries make the rank deficiency frequent
        int64_t range = (iter % 2 == 0) ? 1 : 1000000;
        MyMatrix<int64_t> M_signed(nr, nc);
        ModMatrix M(nr, nc);
        for (int i = 0; i < nr; i++) {
          for (int j = 0; j < nc; j++) {
            M_signed(i, j) = rnd(-range, range);
            M(i, j) = ResModSigned(M_signed(i, j), p);
          }
        }
        if (ModMatrixOfMatrix(M_signed, p) != M) {
          fail("ModMatrixOfMatrix");
        }
        RowSpaceModP rs = ComputeRowSpaceModP(M, p);
        size_t rank = rs.pivots.size();
        if (static_cast<size_t>(rs.basis.rows()) != rank) {
          fail("the number of basis rows");
        }
        // the echelon shape
        for (size_t i = 0; i < rank; i++) {
          for (size_t j = 0; j < rank; j++) {
            if (rs.basis(static_cast<int>(j), static_cast<int>(rs.pivots[i])) !=
                ((i == j) ? 1 : 0)) {
              fail("echelon shape");
            }
          }
        }
        // every row from its coordinates
        for (int i = 0; i < nr; i++) {
          ModVector row(nc);
          for (int j = 0; j < nc; j++) {
            row[j] = M(i, j);
          }
          ModVector c = CoordinatesInRowSpaceModP(rs, row);
          for (int j = 0; j < nc; j++) {
            T_residue s = 0;
            for (size_t u = 0; u < rank; u++) {
              s = ResAdd(s, ResMul(c[u], rs.basis(static_cast<int>(u), j), p),
                         p);
            }
            if (s != row[j]) {
              fail("a row is not recovered from its coordinates");
            }
          }
        }
        // the rank of a square matrix against the determinant
        if (nr == nc && p < (T_small_prime(1) << 31)) {
          std::vector<int64_t> a;
          for (int i = 0; i < nr; i++) {
            for (int j = 0; j < nc; j++) {
              a.push_back(ResToSigned(M(i, j)));
            }
          }
          bool full = DeterminantMatModInt64(a, nr, ResToSigned(p), true) != 0;
          if (full != (static_cast<int>(rank) == nr)) {
            fail("the rank disagrees with the determinant");
          }
        }
        if (RankModPrime(M, p) != rank) {
          fail("RankModPrime");
        }
        // the free kernel basis of the system M x = 0 (the rows of M are
        // the equations): of dimension nc - rank, and M x = 0
        ModMatrix F = FreeKernelBasisModP(rs, nc, p);
        if (static_cast<size_t>(F.rows()) + rank != static_cast<size_t>(nc) ||
            RankModPrime(F, p) != static_cast<size_t>(F.rows())) {
          fail("the dimension of the free kernel basis");
        }
        if (!ModMatrixProduct(M, ModMatrix(F.transpose()), p).isZero()) {
          fail("M x is not zero for the free kernel basis");
        }
        // the left nullspace
        ModMatrix K = LeftNullspaceModP(M, p);
        if (static_cast<size_t>(K.rows()) + rank != static_cast<size_t>(nr) ||
            RankModPrime(K, p) != static_cast<size_t>(K.rows())) {
          fail("the dimension of the left nullspace");
        }
        if (!ModMatrixProduct(K, M, p).isZero()) {
          fail("x M is not zero");
        }
        // x B = v on the echelon basis, v in the row space and maybe not
        if (rank > 0) {
          ModVector x0(rank), v(nc, 0);
          for (auto &x : x0) {
            x = static_cast<T_residue>(rnd(0, ResToSigned(p) - 1));
          }
          for (size_t i = 0; i < rank; i++) {
            for (int j = 0; j < nc; j++) {
              v[j] = ResAdd(v[j],
                            ResMul(x0[i], rs.basis(static_cast<int>(i), j), p),
                            p);
            }
          }
          std::optional<ModVector> opt = SolveRowsModP(rs.basis, v, p);
          if (!opt || *opt != x0) {
            fail("SolveRowsModP on the row space");
          }
          if (static_cast<int>(rank) < nc) {
            ModVector w = v;
            for (int j = 0; j < nc; j++) {
              bool is_pivot = false;
              for (auto piv : rs.pivots) {
                if (static_cast<int>(piv) == j) {
                  is_pivot = true;
                }
              }
              if (!is_pivot) {
                w[j] = ResAdd(w[j], T_residue(1), p);
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
    // the residue operations and the powers
    for (int iter = 0; iter < 2000; iter++) {
      int64_t m_max = (iter % 2 == 0) ? (int64_t(1) << 40) : (int64_t(1) << 62);
      T_residue m = static_cast<T_residue>(rnd(2, m_max));
      int64_t a_s = rnd(-(int64_t(1) << 62), int64_t(1) << 62);
      int64_t b_s = rnd(-(int64_t(1) << 62), int64_t(1) << 62);
      __int128 m128 = static_cast<__int128>(m);
      auto ref_mod = [&](__int128 x) -> T_residue {
        return static_cast<T_residue>(((x % m128) + m128) % m128);
      };
      T_residue a = ResModSigned(a_s, m);
      T_residue b = ResModSigned(b_s, m);
      if (a != ref_mod(a_s) || b != ref_mod(b_s)) {
        fail("ResModSigned");
      }
      if (ResModOf(mpz_class(static_cast<long>(a_s)), m) != a) {
        fail("ResModOf");
      }
      if (ResMul(a, b, m) != ref_mod(static_cast<__int128>(a_s) * b_s)) {
        fail("ResMul");
      }
      if (ResSub(a, b, m) != ref_mod(static_cast<__int128>(a_s) - b_s) ||
          ResAdd(a, b, m) != ref_mod(static_cast<__int128>(a_s) + b_s)) {
        fail("ResAdd or ResSub");
      }
      uint64_t e = static_cast<uint64_t>(rnd(0, 60));
      T_residue naive = 1 % m;
      for (uint64_t i = 0; i < e; i++) {
        naive = ResMul(naive, a, m);
      }
      if (ResPow(a, e, m) != naive) {
        fail("ResPow");
      }
      if (std::gcd(a, m) == 1 && ResMul(a, ResInv(a, m), m) != 1 % m) {
        fail("ResInv");
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
