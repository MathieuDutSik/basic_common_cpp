// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_QUADRATIC_RESIDUE_H_
#define SRC_NUMBER_QUADRATIC_RESIDUE_H_

// clang-format off
#include "TemplateTraits.h"
#include "factorizations.h"
#include <limits>
#include <map>
#include <utility>
#include <vector>
// clang-format on

/*
  References used below:
  OB: Oswald Baumgart, "The Quadratic Reciprocity Law: A collection of classical
  proofs" Edited and translated by Franz Lemmermeyer

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
  See for details: https://en.wikipedia.org/wiki/Jacobi_symbol

  For the Jacobi-symbol, we have (P / Q) = -1 implies that there is no
  quadratic residue.
  But for (P / Q) = 1, we cannot conclude.
 */

/*
  What about computing quadratic residue x^2 = a (mod m) if we know
  the factorization of m.
  If we have m = p1^m1 .... pk^mk
  then we can resolve for each of the prime powers qi = pi^mi.
  Then if for any of those powers there is no solution, then there is
  no solution.
  But if there is a solution for each of them then we can use the
  Chinese remainder theorem to find a X such that
  X = xi (mod qi) for all i.
  Then we have for all i X^2 = a (mod qi) and so X^2 = a (mod m)
  ---
  Now for resolving the equation X^2 = a (mod qi) we can iteratively
  solve mod pi and that will work.
 */

// This is an exhaustive search that works even if m is not prime.
// The function returns a x such that x^2 = a (mod m) if it exists.
template <typename T>
std::optional<T> find_quadratic_residue_exhaustive_kernel(T const &a, T const &m) {
  static_assert(is_implementation_of_Z<T>::value, "Requires T to be a Z ring");
  T two(2);
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
    std::cerr << "QUADRES: x=" << x << " xSqr=" << xSqr << " a=" << a
              << " m=" << m << "\n";
#endif
    if (xSqr == a) {
      return x;
    }
    xSqr += TwoXpOne;
    xSqr = ResInt(xSqr, m);
    TwoXpOne += 2;
    x += 1;
  }
  return {};
}

// Compute the quadratic residue by mapping to a faster numeric if available.
template <typename T, typename Tcomp>
std::optional<T> find_quadratic_residue_exhaustive_Tcomp(T const &a, T const &m) {
  Tcomp a_comp = UniversalScalarConversion<Tcomp, T>(a);
  Tcomp m_comp = UniversalScalarConversion<Tcomp, T>(m);
  std::optional<Tcomp> opt = find_quadratic_residue_exhaustive_kernel(a_comp, m_comp);
  if (opt) {
    Tcomp val_comp = *opt;
    T val = UniversalScalarConversion<T, Tcomp>(val_comp);
    return val;
  } else {
    return {};
  }
}

template <typename T>
std::optional<T> find_quadratic_residue_exhaustive(T const &a_in, T const &m_in) {
  T m = T_abs(m_in);
  T a = ResInt(a_in, m);
  T max16_A = UniversalScalarConversion<T, int16_t>(
      std::numeric_limits<int16_t>::max());
  T max32_A = UniversalScalarConversion<T, int32_t>(
      std::numeric_limits<int32_t>::max());
  T max64_A = UniversalScalarConversion<T, int64_t>(
      std::numeric_limits<int64_t>::max());
  T four(4);
  T max16_B = QuoInt(max16_A, four);
  T max32_B = QuoInt(max32_A, four);
  T max64_B = QuoInt(max64_A, four);
  if (m < max16_B) {
    return find_quadratic_residue_exhaustive_Tcomp<T, int16_t>(a, m);
  }
  if (m < max32_B) {
    return find_quadratic_residue_exhaustive_Tcomp<T, int32_t>(a, m);
  }
  if (m < max64_B) {
    return find_quadratic_residue_exhaustive_Tcomp<T, int64_t>(a, m);
  }
  return find_quadratic_residue_exhaustive_kernel(a, m);
}

template <typename T>
std::optional<T> find_quadratic_residue_map(T const &a_in, std::map<T,size_t> const &m_map) {
  std::vector<T> a;
  std::vector<T> m;
  T m_prod(1);
#ifdef DEBUG_QUADRATIC_RESIDUE
  for (auto & kv: m_map) {
    std::cerr << "QUADRES: p=" << kv.first << " mult=" << kv.second << "\n";
  }
#endif
  for (auto & kv: m_map) {
    T const& p = kv.first;
    T prod = p;
    for (size_t u=1; u<kv.second; u++) {
      prod *= p;
    }
    std::optional<T> opt = find_quadratic_residue_exhaustive(a_in, prod);
    if (opt) {
      T val = *opt;
      a.push_back(val);
      m.push_back(prod);
    } else {
      return {};
    }
    m_prod *= prod;
  }
  T x1 = chinese_remainder_theorem(a, m);
  T x2 = ResInt(x1, m_prod);
#ifdef DEBUG_QUADRATIC_RESIDUE
  T diff = x2 * x2 - a_in;
  T res = ResInt(diff, m_prod);
  if (res != 0) {
    std::cerr << "QUADRES: a=";
    for (auto & val : a) {
      std::cerr << " " << val;
    }
    std::cerr << "\n";
    std::cerr << "QUADRES: m=";
    for (auto & val : m) {
      std::cerr << " " << val;
    }
    std::cerr << "\n";
    std::cerr << "QUADRES: x2 is not a solution\n";
    throw TerminalException{1};
  }
#endif
  return x2;
}



template <typename T>
std::optional<T> find_quadratic_residue(T const &a_in, T const &m_in) {
#ifdef DEBUG_QUADRATIC_RESIDUE
  std::cerr << "QUADRES: a_in=" << a_in << " m_in=" << m_in << "\n";
#endif
  if (T_abs(m_in) == 1) {
    return 0;
  }
  std::map<T, size_t> map = FactorsIntMap(T_abs(m_in));
  return find_quadratic_residue_map(a_in, map);
}



template <typename T>
std::pair<size_t, T> decompose_even_power_odd(T const &val) {
  size_t power = 0;
  T two(2);
  T work_val = val;
  while (true) {
    std::pair<T, T> pair = ResQuoInt(work_val, two);
#ifdef DEBUG_QUADRATIC_RESIDUE
    std::cerr << "work_val=" << work_val << " res=" << pair.first
              << " q=" << pair.second << "\n";
#endif
    if (pair.first > 0) {
      return {power, work_val};
    }
    power += 1;
    work_val = pair.second;
  }
}

template <typename T>
int get_legendre_symbol_power_two(size_t const &power, T const &P) {
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
  T gcd = T_abs(GcdPair(a_in, m));
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
  while (true) {
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
    std::cerr << "Update 1: sign_pq=" << sign_pq << " Pm1d2=" << Pm1d2
              << " Qm1d2=" << Qm1d2 << "\n";
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

template <typename T>
bool is_quadratic_residue_map(T const &a, std::map<T, size_t> const &m) {
  if (m.size() == 0) {
    return true;
  }
  std::optional<T> opt = find_quadratic_residue_map(a, m);
  return opt.has_value();
}

template <typename T>
bool is_quadratic_residue_exhaustive(T const &a, T const &m) {
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
  T gcd = T_abs(GcdPair(a, m));
  if (gcd == 1) {
    return is_quadratic_residue_kernel(a, m);
  }
  return is_quadratic_residue_exhaustive(a, m);
}

// clang-format off
#endif  // SRC_NUMBER_QUADRATIC_RESIDUE_H_
// clang-format on
