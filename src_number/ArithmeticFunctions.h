// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_ARITHMETICFUNCTIONS_H_
#define SRC_NUMBER_ARITHMETICFUNCTIONS_H_

// clang-format off
#include "NumberTheoryCommon.h"
#include "NumberTheoryModular.h"
#include "factorizations.h"
#include <map>
#include <numeric>
#include <optional>
#include <utility>
#include <vector>
// clang-format on

/*
  The arithmetic functions of the integers read off their factorization
  (FactorsIntMap, factorizations.h), for any integer type T:
  - EulerPhi(n), Moebius(n) for n >= 1;
  - IsSquareFree(n), SquareFreePart(n) (with the sign of n) for n != 0;
  - PrimeDivisors(n), OddPrimeDivisors(n): the primes dividing n != 0,
    without multiplicity and sorted;
  - PrimePowerDecomposition(q): the prime p and the exponent k >= 1 with
    q = p^k, nothing when q is not such a power;
  - BinomialCoefficient<T>(n, k), computed in T from integers n, k of
    another type, and the Gaussian binomial [a b]_q.
  The factorization is the expensive part, so a caller needing several of
  them for the same n can factor once and read the map.

  And the small integer functions that need no factorization: the
  valuations with an int64_t prime, the Legendre symbol modulo q, the
  multiplicative order, the Ramanujan sums, the fundamental
  discriminants, the bit length, the integer square root, the next
  prime, and (-1)^e.
 */

template <typename T> T IntegerAbsoluteValue(T const &n) {
  return (n < 0) ? T(-n) : n;
}

template <typename T> T EulerPhi(T const &n) {
  T res(1);
  if (n <= 1) {
    return res;
  }
  for (auto const &[p, mult] : FactorsIntMap(n)) {
    res *= p - 1;
    for (size_t i = 1; i < mult; i++) {
      res *= p;
    }
  }
  return res;
}

template <typename T> int Moebius(T const &n) {
  if (n <= 1) {
    return 1;
  }
  int mu = 1;
  for (auto const &[p, mult] : FactorsIntMap(n)) {
    if (mult > 1) {
      return 0;
    }
    mu = -mu;
  }
  return mu;
}

template <typename T> bool IsSquareFree(T const &n) {
  T n_abs = IntegerAbsoluteValue(n);
  if (n_abs <= 1) {
    return true;
  }
  for (auto const &[p, mult] : FactorsIntMap(n_abs)) {
    if (mult > 1) {
      return false;
    }
  }
  return true;
}

// The square-free part of n != 0, with the sign of n.
template <typename T> T SquareFreePart(T const &n) {
  T res(1);
  T n_abs = IntegerAbsoluteValue(n);
  if (n_abs > 1) {
    for (auto const &[p, mult] : FactorsIntMap(n_abs)) {
      if (mult % 2 == 1) {
        res *= p;
      }
    }
  }
  return (n < 0) ? T(-res) : res;
}

template <typename T> std::vector<T> PrimeDivisors(T const &n) {
  std::vector<T> res;
  T n_abs = IntegerAbsoluteValue(n);
  if (n_abs > 1) {
    for (auto const &[p, mult] : FactorsIntMap(n_abs)) {
      res.push_back(p);
    }
  }
  return res;
}

template <typename T> std::vector<T> OddPrimeDivisors(T const &n) {
  std::vector<T> res;
  for (auto const &p : PrimeDivisors(n)) {
    if (p != 2) {
      res.push_back(p);
    }
  }
  return res;
}

template <typename T>
std::optional<std::pair<T, size_t>> PrimePowerDecomposition(T const &q) {
  if (q <= 1) {
    return {};
  }
  std::map<T, size_t> map = FactorsIntMap(q);
  if (map.size() != 1) {
    return {};
  }
  return *map.begin();
}

// C(n, k) = (n / k) C(n - 1, k - 1), 0 for k < 0.
template <typename T, typename Tint>
T BinomialCoefficient(Tint const &n, Tint const &k) {
  T res(1);
  if (k < 0) {
    return T(0);
  }
  for (Tint i = 1; i <= k; i++) {
    T num = n - (k - i);
    T den = i;
    res *= num;
    res /= den;
  }
  return res;
}

/*
  The Gaussian binomial [a b]_q = prod_{i=1}^{b} (q^(a-b+i) - 1) / (q^i - 1),
  0 unless 0 <= b <= a, by the symmetry [a b] = [a a-b] with b <= a / 2:
  the powers of q are carried along and the quotient is taken once.
 */
template <typename T> T GaussianBinomial(int a, int b, T const &q) {
  if (b < 0 || b > a) {
    return T(0);
  }
  if (2 * b > a) {
    b = a - b;
  }
  T q_num(1);
  for (int i = 0; i < a - b; i++) {
    q_num *= q;
  }
  T q_den(1), num(1), den(1);
  for (int i = 1; i <= b; i++) {
    q_num *= q;
    q_den *= q;
    num *= q_num - 1;
    den *= q_den - 1;
  }
  return num / den;
}

// v_p(x) for a nonzero integer x.
template <typename Tint> int IntegerValuation(Tint const &x, int64_t p) {
  return static_cast<int>(
      prime_valuation(x, UniversalScalarConversion<Tint, int64_t>(p)));
}

// v_p(x) for a nonzero rational x, of either sign.
template <typename T> int RationalValuation(T const &x, int64_t p) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  return rational_prime_valuation(x,
                                  UniversalScalarConversion<Tint, int64_t>(p));
}

// The Legendre symbol (a / q), q an odd prime.
template <typename Tint> int LegendreModQ(Tint const &a, int64_t q) {
  return jacobi_symbol(a, UniversalScalarConversion<Tint, int64_t>(q));
}

// The smallest prime above q.
inline int NextPrimeAfter(int q) {
  int cand = q + 1;
  while (true) {
    bool is_prime = cand >= 2;
    for (int d = 2; d * d <= cand; d++) {
      if (cand % d == 0) {
        is_prime = false;
        break;
      }
    }
    if (is_prime) {
      return cand;
    }
    cand++;
  }
}

// The order of q modulo m, q prime to m.
inline int MultiplicativeOrder(int64_t q, int64_t m) {
  if (m == 1) {
    return 1;
  }
  int64_t x = q % m;
  int o = 1;
  while (x != 1) {
    x = (x * q) % m;
    o++;
  }
  return o;
}

// Whether t is a power of x modulo m (x prime to m; always for m <= 2).
inline bool IsPowerModulo(int64_t t, int64_t x, int64_t m) {
  if (m <= 2) {
    return true;
  }
  int64_t target = ((t % m) + m) % m;
  int64_t y = 1;
  int ord = MultiplicativeOrder(x, m);
  for (int k = 0; k < ord; k++) {
    if (y == target) {
      return true;
    }
    y = (y * (x % m)) % m;
  }
  return false;
}

// The Ramanujan sum c_n(j) = Tr_{Q(zeta_n) / Q}(zeta_n^j) = sum_{e | gcd(n,
// j)} mu(n / e) e.
inline int64_t RamanujanSum(int64_t n, int64_t j) {
  int64_t g = std::gcd(n, ((j % n) + n) % n);
  if (g == 0) {
    g = n;
  }
  int64_t res = 0;
  for (int64_t e = 1; e <= g; e++) {
    if (g % e == 0) {
      res += Moebius(n / e) * e;
    }
  }
  return res;
}

// The discriminant of Q(sqrt(x)) for x non-zero; 1 when x is a square.
inline int64_t FundamentalDiscriminant(int64_t x) {
  int64_t d = SquareFreePart(x);
  if (d == 1) {
    return 1;
  }
  int64_t r = ((d % 4) + 4) % 4;
  return (r == 1) ? d : 4 * d;
}

inline bool IsFundamentalDiscriminant(int64_t D) {
  int64_t r = ((D % 4) + 4) % 4;
  if (r == 1) {
    return IsSquareFree(D);
  }
  if (r == 0) {
    int64_t k = D / 4;
    int64_t rk = ((k % 4) + 4) % 4;
    return (rk == 2 || rk == 3) && IsSquareFree(k);
  }
  return false;
}

// The number of bits of |x|.
template <typename Tint> size_t BitLength(Tint x) {
  if (x < 0) {
    x = -x;
  }
  size_t n = 0;
  Tint two(2);
  while (x > 0) {
    x = QuoInt(x, two);
    n++;
  }
  return n;
}

// The floor of the square root of a nonnegative integer.
template <typename Tint> Tint IntegerSquareRoot(Tint const &a) {
  Tint lo(0), hi(1);
  while (hi * hi <= a) {
    hi *= 2;
  }
  // lo^2 <= a < hi^2
  while (hi - lo > 1) {
    Tint mid = (lo + hi) / 2;
    if (mid * mid <= a) {
      lo = mid;
    } else {
      hi = mid;
    }
  }
  return lo;
}

// (-1)^r x, (-1)^e and the floor of a / 2.
inline int64_t SignedPower(int r, int64_t x) { return (r % 2 == 0) ? x : -x; }
inline int MinusOnePower(int e) { return (e % 2 == 0) ? 1 : -1; }
inline int FloorDiv2(int a) { return (a >= 0) ? a / 2 : -((-a + 1) / 2); }

// clang-format off
#endif  // SRC_NUMBER_ARITHMETICFUNCTIONS_H_
// clang-format on
