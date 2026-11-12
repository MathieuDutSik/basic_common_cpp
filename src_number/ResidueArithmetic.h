// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_GENERA_RESIDUEARITHMETIC_H_
#define SRC_GENERA_RESIDUEARITHMETIC_H_

// clang-format off
#include "NumberTheory.h"
#include "ExceptionsFunc.h"
#include <cstdint>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
// clang-format on

/*
  Residues modulo a word-sized modulus.

  The glue modules, the multiplicity factors and their orbits work with
  residues modulo prime powers q^k. They are kept in machine integers, as
  they are the inner loops of the listings, but the code using them never
  multiplies or adds them directly: every operation below takes residues
  in [0, m) for a modulus 1 <= m <= kMaxResidueModulus and returns one,
  the products going through __int128, so that no operation can overflow.

  The size is checked in one place only, where a modulus is created:
  ResidueModulusOf (from an unbounded integer) and ResiduePrimePower (p^k)
  stop with an error beyond kMaxResidueModulus. With m <= 2^62, a + b < m
  + m <= 2^63 fits, and 2 m (the modulus of a quadratic form on a module of
  exponent m / 2) is a valid modulus when m <= 2^61.
 */

int64_t const kMaxResidueModulus = int64_t(1) << 62;

[[noreturn]] inline void ResidueModulusTooLarge(std::string const &what) {
  std::cerr << "RESIDUE: the modulus " << what << " exceeds 2^62\n";
  throw TerminalException{1};
}

// p^k as a modulus, p >= 2, k >= 0.
inline int64_t ResiduePrimePower(int64_t p, int k) {
  int64_t res = 1;
  for (int i = 0; i < k; i++) {
    if (res > kMaxResidueModulus / p) {
      ResidueModulusTooLarge(std::to_string(p) + "^" + std::to_string(k));
    }
    res *= p;
  }
  return res;
}

// An unbounded positive integer as a modulus.
template <typename Tint> int64_t ResidueModulusOf(Tint const &m) {
  Tint max_T = UniversalScalarConversion<Tint, int64_t>(kMaxResidueModulus);
  if (m < Tint(1)) {
    std::cerr << "RESIDUE: the modulus " << m << " is not positive\n";
    throw TerminalException{1};
  }
  if (m > max_T) {
    std::ostringstream os;
    os << m;
    ResidueModulusTooLarge(os.str());
  }
  return UniversalScalarConversion<int64_t, Tint>(m);
}

// a mod m in [0, m), for any a.
inline int64_t ResMod(int64_t a, int64_t m) {
  int64_t r = a % m;
  return (r < 0) ? r + m : r;
}

// x mod m in [0, m), for an unbounded integer x.
template <typename Tint> int64_t ResModOf(Tint const &x, int64_t m) {
  Tint m_T = UniversalScalarConversion<Tint, int64_t>(m);
  Tint r = ResInt(x, m_T);
  if (r < Tint(0)) {
    r += m_T;
  }
  return UniversalScalarConversion<int64_t, Tint>(r);
}

// The operations on residues a, b in [0, m).
inline int64_t ResAdd(int64_t a, int64_t b, int64_t m) {
  int64_t s = a + b;
  return (s >= m) ? s - m : s;
}

inline int64_t ResSub(int64_t a, int64_t b, int64_t m) {
  int64_t s = a - b;
  return (s < 0) ? s + m : s;
}

inline int64_t ResNeg(int64_t a, int64_t m) { return (a == 0) ? 0 : m - a; }

inline int64_t ResMul(int64_t a, int64_t b, int64_t m) {
  __int128 prod = static_cast<__int128>(a) * static_cast<__int128>(b);
  return static_cast<int64_t>(prod % static_cast<__int128>(m));
}

// The inverse of a modulo m, an error when a is not invertible. The
// coefficients of the extended Euclid algorithm are bounded by m.
inline int64_t ResInv(int64_t a, int64_t m) {
  int64_t r0 = m, r1 = ResMod(a, m);
  int64_t s0 = 0, s1 = 1;
  while (r1 != 0) {
    int64_t q = r0 / r1;
    int64_t r2 = r0 - q * r1;
    r0 = r1;
    r1 = r2;
    int64_t s2 = s0 - q * s1;
    s0 = s1;
    s1 = s2;
  }
  if (r0 != 1) {
    std::cerr << "RESIDUE: " << a << " is not invertible modulo " << m << "\n";
    throw TerminalException{1};
  }
  return ResMod(s0, m);
}

// The product of sizes, saturated at the largest size_t, for the bounds of
// the listings: a product too large to be listed is never mistaken for a
// small one.
inline size_t SaturatedSizeProduct(size_t a, size_t b) {
  if (a != 0 && b > std::numeric_limits<size_t>::max() / a) {
    return std::numeric_limits<size_t>::max();
  }
  return a * b;
}

// clang-format off
#endif  // SRC_GENERA_RESIDUEARITHMETIC_H_
// clang-format on
