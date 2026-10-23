// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_MatrixMod.h"
#include "MAT_MatrixModInt64.h"
#include <string>
#include <vector>
// clang-format on

/*
  Tests of the small dense matrices modulo m against independent
  computations: the determinant against the cofactor expansion (any m) and
  against DeterminantMatMod (prime m), the inverse by A A^{-1} = I, and its
  existence exactly when the determinant is a unit modulo m.
 */

int64_t det_cofactor(std::vector<int64_t> const &a, int n, int64_t m) {
  if (n == 1) {
    return ResidueModPrime(a[0], m);
  }
  int64_t det = 0;
  std::vector<int64_t> minor((n - 1) * (n - 1));
  for (int j = 0; j < n; j++) {
    int pos = 0;
    for (int i2 = 1; i2 < n; i2++) {
      for (int j2 = 0; j2 < n; j2++) {
        if (j2 != j) {
          minor[pos] = a[i2 * n + j2];
          pos += 1;
        }
      }
    }
    int64_t term = (a[j] * det_cofactor(minor, n - 1, m)) % m;
    det = ResidueModPrime((j % 2 == 0) ? det + term : det - term, m);
  }
  return det;
}

bool is_prime(int64_t m) {
  if (m < 2) {
    return false;
  }
  for (int64_t d = 2; d * d <= m; d++) {
    if (m % d == 0) {
      return false;
    }
  }
  return true;
}

int main() {
  HumanTime time;
  try {
    int n_error = 0;
    std::vector<int64_t> l_mod{2, 3, 5, 7, 11, 13, 101, 2147483647, 1, 4, 6, 8, 9, 12, 15, 30, 121, 1000};
    size_t n_invertible = 0, n_singular = 0;
    for (auto &m : l_mod) {
      bool prime = is_prime(m);
      for (int iter = 0; iter < 300; iter++) {
        int n = random_int(1, 6);
        std::vector<int64_t> a(n * n);
        // Small entries make the singular matrices frequent.
        int range = (iter % 2 == 0) ? 2 : 1000;
        for (auto &x : a) {
          x = ResidueModPrime(random_int(-range, range), m);
        }
        int64_t det = DeterminantMatModInt64(a, n, m, prime);
        int64_t det_ref = (m == 1) ? 0 : det_cofactor(a, n, m);
        if (det != det_ref) {
          std::cerr << "m=" << m << " n=" << n << " det=" << det << " det_ref=" << det_ref << "\n";
          n_error += 1;
        }
        if (prime && m < 1000) {
          MyMatrix<mpz_class> A(n, n);
          for (int i = 0; i < n * n; i++) {
            A(i / n, i % n) = mpz_class(static_cast<long>(a[i]));
          }
          mpz_class det_mpz = DeterminantMatMod(A, mpz_class(static_cast<long>(m)));
          if (det_mpz != mpz_class(static_cast<long>(det))) {
            std::cerr << "m=" << m << " det differs from DeterminantMatMod\n";
            n_error += 1;
          }
        }
        std::optional<std::vector<int64_t>> opt = InverseMatModInt64(a, n, m, prime);
        bool unit = (m == 1) || std::gcd(det_ref, m) == 1;
        if (opt.has_value() != unit) {
          std::cerr << "m=" << m << " n=" << n << " the inverse exists=" << opt.has_value()
                    << " but det=" << det_ref << "\n";
          n_error += 1;
          continue;
        }
        if (!opt) {
          n_singular += 1;
          continue;
        }
        n_invertible += 1;
        std::vector<int64_t> const &inv = *opt;
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++) {
            int64_t sum = 0;
            for (int k = 0; k < n; k++) {
              sum = (sum + a[i * n + k] * inv[k * n + j]) % m;
            }
            int64_t expected = (i == j) ? ResidueModPrime(1, m) : 0;
            if (sum != expected) {
              std::cerr << "m=" << m << " n=" << n << " A A^{-1} is not the identity\n";
              n_error += 1;
              i = n;
              break;
            }
          }
        }
      }
    }
    std::cerr << "n_invertible=" << n_invertible << " n_singular=" << n_singular << "\n";
    if (n_error > 0) {
      std::cerr << "Test_MatrixModInt64: " << n_error << " errors\n";
      return 1;
    }
    std::cerr << "Test_MatrixModInt64: all tests passed\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_MatrixModInt64\n";
    exit(e.eVal);
  }
  runtime(time);
}
