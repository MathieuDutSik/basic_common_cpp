// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheory.h"
#ifdef ENABLE_FLINT_SUPPORT
#include "NumberTheoryFlint.h"
#endif
#include "TypeConversion.h"
#include "NumberTheoryModular.h"
#include "PolynomialRealRoots.h"
// clang-format on

/*
  Checks of NumberTheoryModular.h and PolynomialRealRoots.h over the integer
  and rational types: every value is compared with a brute force reference
  computed with int64_t.
 */

int64_t ref_power_modulo(int64_t a, int64_t e, int64_t m) {
  int64_t r = 1 % m;
  int64_t b = ((a % m) + m) % m;
  for (int64_t i = 0; i < e; i++) {
    r = (r * b) % m;
  }
  return r;
}

// the Legendre symbol by Euler's criterion, n an odd prime
int64_t ref_legendre(int64_t a, int64_t p) {
  int64_t r = ref_power_modulo(a, (p - 1) / 2, p);
  if (r == 0) {
    return 0;
  }
  return (r == 1) ? 1 : -1;
}

// the Kronecker symbol from the factorization of n > 0
int64_t ref_kronecker(int64_t a, int64_t n) {
  int64_t res = 1;
  int64_t m = n;
  for (int64_t p = 2; p <= m; p++) {
    while (m % p == 0) {
      m /= p;
      if (p == 2) {
        int64_t r = ((a % 8) + 8) % 8;
        if (r % 2 == 0) {
          return 0;
        }
        if (r == 3 || r == 5) {
          res = -res;
        }
      } else {
        res *= ref_legendre(a, p);
      }
    }
  }
  return res;
}

template <typename T> size_t check_integer_type(std::string const &name) {
  size_t n_error = 0;
  auto conv = [](int64_t x) -> T { return UniversalScalarConversion<T, int64_t>(x); };
  auto back = [](T const &x) -> int64_t { return UniversalScalarConversion<int64_t, T>(x); };
  for (int64_t m = 1; m < 40; m++) {
    for (int64_t a = -30; a < 30; a++) {
      for (int64_t e = 0; e < 12; e++) {
        int64_t v = back(power_modulo(conv(a), conv(e), conv(m)));
        if (v != ref_power_modulo(a, e, m)) {
          std::cerr << name << ": power_modulo(" << a << ", " << e << ", " << m
                    << ") = " << v << "\n";
          n_error++;
        }
      }
    }
  }
  for (int64_t n = 1; n < 120; n++) {
    for (int64_t a = -60; a < 60; a++) {
      int k = kronecker_symbol(conv(a), conv(n));
      if (k != ref_kronecker(a, n)) {
        std::cerr << name << ": kronecker(" << a << ", " << n << ") = " << k
                  << " ref=" << ref_kronecker(a, n) << "\n";
        n_error++;
      }
      if (n % 2 == 1) {
        int j = jacobi_symbol(conv(a), conv(n));
        if (j != ref_kronecker(a, n)) {
          std::cerr << name << ": jacobi(" << a << ", " << n << ") = " << j << "\n";
          n_error++;
        }
      }
    }
  }
  for (int64_t p : {2, 3, 5, 7}) {
    for (int64_t x = -500; x < 500; x++) {
      if (x == 0) {
        continue;
      }
      size_t ref = 0;
      int64_t y = x;
      while (y % p == 0) {
        y /= p;
        ref++;
      }
      if (prime_valuation(conv(x), conv(p)) != ref) {
        std::cerr << name << ": prime_valuation(" << x << ", " << p << ")\n";
        n_error++;
      }
    }
  }
  // exact roots: the y with y^r = x, y >= 0 for r even, by brute force
  for (size_t r = 1; r < 6; r++) {
    for (int64_t x = -300; x <= 300; x++) {
      std::optional<int64_t> ref;
      for (int64_t y = -300; y <= 300; y++) {
        int64_t z = 1;
        for (size_t i = 0; i < r; i++) {
          z *= y;
        }
        if (z == x && (r % 2 == 1 || y >= 0)) {
          ref = y;
        }
      }
      std::optional<T> root = exact_integer_root(conv(x), r);
      bool ok = (root.has_value() == ref.has_value()) &&
                (!ref || back(*root) == *ref);
      if (!ok) {
        std::cerr << name << ": exact_integer_root(" << x << ", " << r << ")\n";
        n_error++;
      }
    }
  }
  std::cerr << name << ": " << n_error << " errors\n";
  return n_error;
}

template <typename Tq, typename Tz> size_t check_rational_type(std::string const &name) {
  size_t n_error = 0;
  auto conv = [](int64_t x) -> Tq { return UniversalScalarConversion<Tq, int64_t>(x); };
  // valuations of rationals
  for (int64_t num = -40; num <= 40; num++) {
    for (int64_t den = 1; den <= 40; den++) {
      if (num == 0) {
        continue;
      }
      Tq x = conv(num);
      x /= conv(den);
      for (int64_t p : {2, 3, 5}) {
        auto vp = [&](int64_t y) -> int {
          int v = 0;
          while (y % p == 0) {
            y /= p;
            v++;
          }
          return v;
        };
        int ref = vp(num) - vp(den);
        int v = rational_prime_valuation(x, UniversalScalarConversion<Tz, int64_t>(p));
        if (v != ref) {
          std::cerr << name << ": rational_prime_valuation(" << num << "/" << den
                    << ", " << p << ") = " << v << " ref=" << ref << "\n";
          n_error++;
        }
      }
    }
  }
  // real roots
  Tq eps = Tq(1);
  for (int i = 0; i < 100; i++) {
    eps /= 2;
  }
  struct Case {
    std::vector<int64_t> poly;
    size_t n_root;
  };
  std::vector<Case> cases = {{{-2, 0, 1}, 2},
                             {{-6, 11, -6, 1}, 3},
                             {{-1, -2, 1, 1}, 3},
                             {{1, 0, -10, 0, 1}, 4},
                             {{1, 0, 1}, 0},
                             {{2, -8, 0, 0, 1}, 2},
                             {{-1, 1}, 1}};
  for (auto const &c : cases) {
    std::vector<Tq> poly;
    for (auto x : c.poly) {
      poly.push_back(conv(x));
    }
    std::vector<Tq> roots = real_roots_approximation(poly, eps);
    if (roots.size() != c.n_root) {
      std::cerr << name << ": " << roots.size() << " roots instead of " << c.n_root << "\n";
      n_error++;
      continue;
    }
    for (auto const &r : roots) {
      Tq rm = r - eps, rp = r + eps;
      Tq lo = polynomial_evaluation(poly, rm);
      Tq hi = polynomial_evaluation(poly, rp);
      Tq prod = lo * hi;
      if (polynomial_evaluation(poly, r) != 0 && prod > 0) {
        std::cerr << name << ": no sign change around a root\n";
        n_error++;
      }
    }
  }
  for (int64_t x = 0; x < 50; x++) {
    Tq s = square_root_approximation(conv(x), eps);
    Tq lo = s - eps, hi = s + eps;
    Tq lo2 = lo * lo, hi2 = hi * hi;
    if (lo > 0 && lo2 > conv(x)) {
      n_error++;
    }
    if (hi2 < conv(x)) {
      n_error++;
    }
  }
  std::cerr << name << ": " << n_error << " errors\n";
  return n_error;
}

int main() {
  try {
    size_t n_error = 0;
    n_error += check_integer_type<int64_t>("int64_t");
    n_error += check_integer_type<mpz_class>("mpz_class");
    n_error += check_integer_type<boost::multiprecision::cpp_int>("cpp_int");
    n_error += check_integer_type<boost::multiprecision::mpz_int>("mpz_int");
#ifdef ENABLE_FLINT_SUPPORT
    n_error += check_integer_type<fmpz_class>("fmpz_class");
#endif
    n_error += check_rational_type<mpq_class, mpz_class>("mpq_class");
    n_error += check_rational_type<boost::multiprecision::cpp_rational,
                                   boost::multiprecision::cpp_int>("cpp_rational");
    n_error += check_rational_type<boost::multiprecision::mpq_rational,
                                   boost::multiprecision::mpz_int>("mpq_rational");
#ifdef ENABLE_FLINT_SUPPORT
    n_error += check_rational_type<fmpq_class, fmpz_class>("fmpq_class");
#endif
    if (n_error > 0) {
      std::cerr << "Test_ModularRoots: " << n_error << " errors\n";
      return 1;
    }
    std::cerr << "Test_ModularRoots: all checks passed\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
