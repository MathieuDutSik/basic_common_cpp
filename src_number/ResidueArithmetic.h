// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_RESIDUEARITHMETIC_H_
#define SRC_NUMBER_RESIDUEARITHMETIC_H_

// clang-format off
#include "NumberTheory.h"
#include "ExceptionsFunc.h"
#include <cstdint>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
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

// A small prime, and a residue modulo a word-sized modulus (a prime or a
// prime power), in [0, m). Both are unsigned: a signed value enters only
// through ResModSigned or ResModOf, and the operations below refuse signed
// arguments (deleted templates), so that a negative number is never
// silently read as a huge residue.
using T_small_prime = uint64_t;
using T_residue = uint64_t;

T_residue const kMaxResidueModulus = T_residue(1) << 62;

[[noreturn]] inline void ResidueModulusTooLarge(std::string const &what) {
  std::cerr << "RESIDUE: the modulus " << what << " exceeds 2^62\n";
  throw TerminalException{1};
}

// p^k as a modulus, p >= 2, k >= 0.
inline T_residue ResiduePrimePower(T_small_prime p, int k) {
  T_residue res = 1;
  for (int i = 0; i < k; i++) {
    if (res > kMaxResidueModulus / p) {
      ResidueModulusTooLarge(std::to_string(p) + "^" + std::to_string(k));
    }
    res *= p;
  }
  return res;
}

// An unbounded positive integer as a modulus.
template <typename Tint> T_residue ResidueModulusOf(Tint const &m) {
  Tint max_T = UniversalScalarConversion<Tint, uint64_t>(kMaxResidueModulus);
  if (m < Tint(1)) {
    std::cerr << "RESIDUE: the modulus " << m << " is not positive\n";
    throw TerminalException{1};
  }
  if (m > max_T) {
    std::ostringstream os;
    os << m;
    ResidueModulusTooLarge(os.str());
  }
  return UniversalScalarConversion<uint64_t, Tint>(m);
}

// a mod m in [0, m) for a signed machine integer a.
inline T_residue ResModSigned(int64_t a, T_residue m) {
  int64_t m_s = static_cast<int64_t>(m);
  int64_t r = a % m_s;
  return static_cast<T_residue>((r < 0) ? r + m_s : r);
}

// x mod m in [0, m), for an unbounded integer x of either sign.
template <typename Tint> T_residue ResModOf(Tint const &x, T_residue m) {
  Tint m_T = UniversalScalarConversion<Tint, uint64_t>(m);
  Tint r = ResInt(x, m_T);
  if (r < Tint(0)) {
    r += m_T;
  }
  return UniversalScalarConversion<uint64_t, Tint>(r);
}

// a residue as a signed machine integer, in [0, m), m <= 2^62
inline int64_t ResToSigned(T_residue a) { return static_cast<int64_t>(a); }

// The operations on residues a, b in [0, m).
inline T_residue ResAdd(T_residue a, T_residue b, T_residue m) {
  T_residue s = a + b;
  return (s >= m) ? s - m : s;
}

inline T_residue ResSub(T_residue a, T_residue b, T_residue m) {
  return (a >= b) ? a - b : a + (m - b);
}

inline T_residue ResNeg(T_residue a, T_residue m) {
  return (a == 0) ? 0 : m - a;
}

inline T_residue ResMul(T_residue a, T_residue b, T_residue m) {
  unsigned __int128 prod =
      static_cast<unsigned __int128>(a) * static_cast<unsigned __int128>(b);
  return static_cast<T_residue>(prod % m);
}

// a^e mod m, by binary exponentiation.
inline T_residue ResPow(T_residue a, uint64_t e, T_residue m) {
  T_residue res = (m == 1) ? 0 : 1;
  while (e > 0) {
    if (e & 1) {
      res = ResMul(res, a, m);
    }
    a = ResMul(a, a, m);
    e >>= 1;
  }
  return res;
}

// The inverse of a modulo m, an error when a is not invertible. The
// coefficients of the extended Euclid algorithm are bounded by m <= 2^62.
inline T_residue ResInv(T_residue a, T_residue m) {
  int64_t r0 = static_cast<int64_t>(m), r1 = static_cast<int64_t>(a % m);
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
  return ResModSigned(s0, m);
}

// No signed argument, not even a literal: a negative value would be read
// as a huge residue (ResMul(a, -1, p) is not -a for an odd p). Reduce it
// with ResModSigned first, and write the constants as T_residue(c).
template <typename T>
concept ResidueSignedArgument = std::is_signed_v<T>;
template <typename A, typename B, typename M>
  requires(ResidueSignedArgument<A> || ResidueSignedArgument<B> ||
           ResidueSignedArgument<M>)
T_residue ResAdd(A, B, M) = delete;
template <typename A, typename B, typename M>
  requires(ResidueSignedArgument<A> || ResidueSignedArgument<B> ||
           ResidueSignedArgument<M>)
T_residue ResSub(A, B, M) = delete;
template <typename A, typename B, typename M>
  requires(ResidueSignedArgument<A> || ResidueSignedArgument<B> ||
           ResidueSignedArgument<M>)
T_residue ResMul(A, B, M) = delete;
template <typename A, typename M>
  requires(ResidueSignedArgument<A> || ResidueSignedArgument<M>)
T_residue ResNeg(A, M) = delete;
template <typename A, typename M>
  requires(ResidueSignedArgument<A> || ResidueSignedArgument<M>)
T_residue ResInv(A, M) = delete;
template <typename A, typename E, typename M>
  requires(ResidueSignedArgument<A> || ResidueSignedArgument<M>)
T_residue ResPow(A, E, M) = delete;

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
#endif  // SRC_NUMBER_RESIDUEARITHMETIC_H_
// clang-format on
