// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_QUADRATIC_RESIDUE_H_
#define SRC_NUMBER_QUADRATIC_RESIDUE_H_

// clang-format off
#include "TemplateTraits.h"
#include <map>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_QUADRATIC_RESIDUE
#endif




/*
  References used below:
  OB: Oswald Baumgart, "The Quadratic Reciprocity Law: A collection of classical proofs"
     Edited and translated by Franz Lemmermeyer

  Th (OB[End of Chapter 2]): If P and Q are positive coprime integers
      then we have
      (P / Q)  (Q / P) = (-1)^( (P-1)/2 . (Q-1)/2 )

  Could we use the theorem to compute the residue in practice?
  * If the numbers are not coprime.
    We want to test the existence of a solution for
    x^2 = P + alpha Q
    But P and Q have a common factor say h.
    That reduces to
    If H has a common square, then we reduce on both sides and all is well.
    BUT: If there are primes like h that occur only one time.
    Then x has to to be divisible by h. That reduces the equation to
    h x^2 = P + alpha Q     with P and Q not divisible by h
    What to do next?
    h is invertible in Z/QZ since it is not a factor.
    Take the invert of h and the equation reduces to
    x^2 = h^{-1} P mod Q
    So, we have indeed a reduction of complexity.
    The number h^{-1} P is coprime with Q since if it were not we would
    not have been able to invert.
  * Now, if 2 occurs as a prime factor of both, we can apply the same
    trick above.
  * We have the fomula
    (2 / P) = (-1)^{ (P^2 - 1)/8 }  See Suppelementary Laws in OB
  * We have (p / 2) = 1 for the odd case.

  But maybe the best is to compute with Kronecker algorithm
  Problem is conditions exposed in 

  
 */





// This is an exhaustive search that works even if m is not prime.
template <typename T>
std::optional<T> find_quadratic_residue(T const &a, T const &m_in) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  T m = T_abs(m_in);
  T two(2);
  T a_mod = ResInt(a, m);
  T res = ResInt(m, two);
  T upper(0);
  if (res == 0) {
    upper = QuoInt(m, two) + 1;
  } else {
    T mP1 = m + 1;
    upper = QuoInt(mP1, two);
  }
  T x(0);
  T TwoXpOne(1);
  T xSqr(0);
#ifdef DEBUG_QUADRATIC_RESIDUE
  std::cerr << "QUADRES: upper=" << upper << "\n";
#endif
  while (x != upper) {
#ifdef DEBUG_QUADRATIC_RESIDUE
    std::cerr << "QUADRES: x=" << x << " xSqr=" << xSqr << " a_mod=" << a_mod
              << " m=" << m << "\n";
#endif
    if (xSqr == a_mod) {
      return x;
    }
    xSqr += TwoXpOne;
    xSqr = ResInt(xSqr, m);
    TwoXpOne += 2;
    x += 1;
  }
  return {};
}

template <typename T> bool is_quadratic_residue(T const &a, T const &m) {
  std::optional<T> opt = find_quadratic_residue(a, m);
  return opt.has_value();
}

// clang-format off
#endif  // SRC_NUMBER_QUADRATIC_RESIDUE_H_
// clang-format on
