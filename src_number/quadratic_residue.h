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

  Combined together that gets us an algorithm.

  But what is computed is not the Legendre symbol but the Jacobi symbol.

  For the Jacobi-symbol, we have (P / Q) = -1 implies that there is no
  quadratic residue.
  But for (P / Q) = 1, we cannot conclude.
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
#ifdef DEBUG_QUADRATIC_RESIDUE_DISABLE
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

template<typename T>
std::pair<size_t, T> decompose_even_power_odd(T const& val) {
  size_t power = 0;
  T two(2);
  T work_val = val;
  while(true) {
    std::pair<T, T> pair = ResQuoInt(work_val, two);
#ifdef DEBUG_QUADRATIC_RESIDUE
    std::cerr << "work_val=" << work_val << " res=" << pair.first << " q=" << pair.second << "\n";
#endif
    if (pair.first > 0) {
      return {power, work_val};
    }
    power += 1;
    work_val = pair.second;
  }
}

template<typename T>
int get_legendre_symbol_power_two(size_t const& power, T const& P) {
  T val = P * P - 1;
  T eight(8);
  T quot = QuoInt(val, eight);
#ifdef DEBUG_QUADRATIC_RESIDUE
  T res = ResInt(val, eight);
  if (res != 0) {
    std::cerr << "The residue is not what we expected\n";
    throw TerminalException{1};
  }
#endif
  T two(2);
  T res2 = ResInt(quot, two);
  if (res2 > 0) {
    size_t power_res = power % 2;
    if (power_res > 0) {
      return -1;
    } else {
      return 1;
    }
  } else {
    return 1;
  }
}

template <typename T> bool compute_jacobi_symbol(T const &a_in, T const &m) {
#ifdef DEBUG_QUADRATIC_RESIDUE
  T gcd = T_abs(GenericGcd(a_in, m));
  if (gcd != 1) {
    std::cerr << "The algorithm is not working for gcd > 1 right now\n";
    throw TerminalException{1};
  }
  if (m < 0) {
    std::cerr << "We need a and b greater than 0\n";
    throw TerminalException{1};
  }
#endif
  T a = ResInt(a_in, m);
  if (a == 0) {
    return true;
  }
  std::pair<size_t, T> pair_a = decompose_even_power_odd(a);
  std::pair<size_t, T> pair_m = decompose_even_power_odd(m);
  if (pair_a.first > 0 && pair_m.first > 0) {
    std::cerr << "We should not reach that problem. Illogic\n";
    throw TerminalException{1};
  }
  T two(2);
  int symbol = get_legendre_symbol_power_two(pair_a.first, pair_m.second);
#ifdef DEBUG_QUADRATIC_RESIDUE
  std::cerr << "power=" << pair_a.first << " P=" << pair_m.second << "\n";
  std::cerr << "symbol=" << symbol << "\n";
#endif
  T P = pair_a.second;
  T Q = pair_m.second;
  while(true) {
    // 1: termination test
    if (P == 1) {
      break;
    }
    // 2: Computing the term (-1)^( (P-1)/2 . (Q-1)/2 )
    T Pm1 = P - 1;
    T Qm1 = Q - 1;
    T Pm1d2 = QuoInt(Pm1, two);
    T Qm1d2 = QuoInt(Qm1, two);
    T res2_P = ResInt(Pm1d2, two);
    T res2_Q = ResInt(Qm1d2, two);
    int sign_pq = 1;
    if (res2_P > 0 && res2_Q > 0) {
      sign_pq = -1;
    }
    symbol *= sign_pq;
#ifdef DEBUG_QUADRATIC_RESIDUE
    std::cerr << "Update 1: sign_pq=" << sign_pq << " Pm1d2=" << Pm1d2 << " Qm1d2=" << Qm1d2 << "\n";
#endif
    // 3: Computing the residue
    T Qres = ResInt(Q, P);
    Q = P;
    // 4: Now decomposing the power of two
    std::pair<size_t, T> pair = decompose_even_power_odd(Qres);
    // 5: Computing the term from the power of two
    int sign_two = get_legendre_symbol_power_two(pair.first, P);
#ifdef DEBUG_QUADRATIC_RESIDUE
    std::cerr << "Update 1: sign_two=" << sign_two << "\n";
#endif
    // 6: Updating
    symbol *= sign_two;
    P = pair.second;
#ifdef DEBUG_QUADRATIC_RESIDUE
    std::cerr << "Now P=" << P << " Q=" << Q << "\n";
#endif
  }
  if (symbol == 1) {
    return true;
  } else {
    return false;
  }
}





template <typename T> bool is_quadratic_residue_exhaustive(T const &a, T const &m) {
  std::optional<T> opt = find_quadratic_residue(a, m);
  return opt.has_value();
}

template <typename T> bool is_quadratic_residue_kernel(T const &a, T const &m) {
  bool test_jacobi = compute_jacobi_symbol(a, m);
#ifdef DEBUG_QUADRATIC_RESIDUE
  bool test_exhaust = is_quadratic_residue_exhaustive(a, m);
  if (test_exhaust && !test_jacobi) {
    std::cerr << "test_exhaust=" << test_exhaust << "\n";
    std::cerr << "test_jacobi=" << test_jacobi << "\n";
    std::cerr << "That is not what we expect from the Jacobi symbol\n";
    std::cerr << "incoherency in the result\n";
    throw TerminalException{1};
  }
  return test_exhaust;
#else
  if (!test_jacobi) {
    return false;
  }
  return is_quadratic_residue_exhaustive(a, m);
#endif
}

template <typename T> bool is_quadratic_residue(T const &a, T const &m) {
  T gcd = T_abs(GenericGcd(a, m));
  if (gcd == 1) {
    return is_quadratic_residue_kernel(a, m);
  }
  return is_quadratic_residue_exhaustive(a, m);
}

// clang-format off
#endif  // SRC_NUMBER_QUADRATIC_RESIDUE_H_
// clang-format on
