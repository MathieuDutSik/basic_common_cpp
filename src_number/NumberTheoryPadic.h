// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYPADIC_H_
#define SRC_NUMBER_NUMBERTHEORYPADIC_H_

// clang-format off
#include "Temp_common.h"
#include "InputOutput.h"
#include "quadratic_residue.h"
#include <limits>
#include <string>
#include <algorithm>
#include <utility>
#include <vector>
// clang-format on

#ifdef DEBUG
#define DEBUG_PADIC
#endif

/*
  Computation with p-adic numbers.

  ---

  Computing with p-adic numbers.
  We have the degree and a vector of entries of length d with all coefficients
  between 0 and p-1. Then:
  * Addition looks trivial, but the question of relative exponent
    muddies things considerably.
  * Multiplication is two loops though we have to care about the
  * Inverse is done in the following way:
    - The exponent issue is trivial.
    - We write x as (x_0, x_1, ...., x_d)
    - We write the inverse y = x^{-1} as (y_0, y_1, ...., y_d)
      The y_0 is found by inverting modulo p and lifting back to Z.
    - Then the y_1 is found via
      1 = x y
        = (x_0 + p x_1) (y_0 + p y_1)
        = x_0 y_0 + p (x_0 y_1 + x_1 y_0)
      Expressing x_0 y_0 = 1 + p h and expanding we can effectively find the
  inverse.
    - There is a more straightforward method. Expressing the two terms
      together we get
      xy = 1 + p^d h
      So, by using the Gcd, we can effectively do the computation.

  ---

  We have R_p = Q_p^* / (Q_p^*)^2
  This is for interest in our computation:
  * For p>2 we have R_p isomorphic to Z_2 x Z_2.
    The two components can be seen from:
    - The square is multiplying the valuation by two. So, elements of odd
      valuations cannot be square.
    - Then we have the square root from the quadratic residue.
    Then we apply the Hensel lemma for getting it.
  * We have R_2 = Z_2 x Z_2 x Z_2
    The components are coming in the following way:
    - The square is multiplying the valuation by two. So, elements of odd
      valuations cannot be square.
    - The square of x = 1 + 2 u is
      x^2 = 1 + 4u + 4u^2
          = 1 + 4 (u + u^2)
      So, the coefficient of 2 is zero, which gives another restriction.
    - Then the entry u + u^2 is always even, so that gives another factor 2.
    Therefore, if we admit the result that R_2 has size 8, then we have an
  algorithm for checking See
  https://math.stackexchange.com/questions/473595/characterization-of-integers-which-has-a-2-adic-square-root
  So, if we can compute in p-adic up to some level, then we can decide whether
  entries are trivial in R_p or not.

  ---

  We do not try to git the P-adic with other data types since a constructor
  Padic x(1);
  would make sense if we were always working with the same p. In that case the
  "p" would work the same way as the d in the quadratic field Q(sqrt(d)).

  But that is precisely not our use case. So, instead we compute with the p put
  as argument.

  We encode a P-adic number by the following data structure
  struct Padic {
    int eff_valuation;
    size_t accuracy;
    std::vector<T> coefficients;
  }
  So, the corresponding number will be
  x = sum_{i=0}^{coefficients.size()-1} coefficients[i] p^{i + eff_valuation}

  Thus a zero has an effective valuation of zero, an accuracy of infinity
  From that type, we can build the functions:
  * valuation() that will be always greater than eff_valuation of possibly
  returns std::numeric_limits<size_t>::max() if zero. If there is not enough
  coefficients to conclude then throw a PadicPrecisionException{1};
  * coeff(i) returns the coefficient. Throw a PadicPrecisionException{2} if not
  accessible.

  Now about implementation:
  * We need a function accuracy_reduction that reduce the number of coefficients
  and the accuracy if non-zero. For zero entries (ones with accuracy = MAX), the
  accuracy is not tested.
  * We need a function is_zero that tests if a number is zero or not.
  * The computation of product returns object whose accuracy is the minimum of
  both accuracies the eff_valuation is computed as the sum of both accuracies.
  * The inverse is done in the following way:
    - Same accuracy.
    - Assumes that the first coefficient is non-zero (that is eff_valuation =
  valuation). If zero, then throw a PadicPrecisionException{3}
    - in return eff_valuation = - eff_valuation
  * The sum is the most problematic.
    - We need to take the eff_valuation and compute what can be computed.

  We are not trying to be partically efficient, just something reasonable for
  computing. The PadicPrecisionException allows for increase of precision if the
  user can solve his problem.
 */

// Running with P-adic number is always done for a precision that is finite
// which is intrinsically a problem because sometimes that just does not work.
// Possible values:
// 1: Valuation exception. Cannot find a non-zero entry, so cannot get the
// valuation. 2: Coefficient error. Cannot access the coefficient that is
// sought. 3: Inverse error. The first coefficient has to be non-zero. 4:
// precision error. Trying to compute inverse with infinite precision, which is
// impossible on a computer. 5: square root test. For p=2 we need 3 coefficients
// to compute. Otherwise, just one number suffices.

struct PadicPrecisionException {
  int val;
};

template <typename T> struct Padic {
  int eff_valuation;
  size_t precision;
  std::vector<T> coefficients;
};

template <typename T> std::string Padic_to_string(Padic<T> const &x) {
  std::stringstream s;
  s << "eff_valuation=" << x.eff_valuation << " precision=" << x.precision
    << " coefficients=";
  for (auto &val : x.coefficients) {
    s << " " << val;
  }
  std::string converted(s.str());
  return converted;
}

template <typename T>
void Padic_debug_print(Padic<T> const &x, std::ostream &os) {
  os << Padic_to_string(x) << "\n";
}

template <typename T>
Padic<T> Padic_from_positive_integer(T const &val, T const &p) {
  std::vector<T> coefficients;
  T val_work = val;
  int eff_valuation = 0;
  bool has_nz = false;
  while (true) {
    if (val_work == 0) {
      break;
    }
    std::pair<T, T> pair = ResQuoInt(val_work, p);
    if (pair.first == 0) {
      if (!has_nz) {
        eff_valuation += 1;
      } else {
        coefficients.push_back(pair.first);
      }
    } else {
      has_nz = true;
      coefficients.push_back(pair.first);
    }
    val_work = pair.second;
  }
  // We have a full integer as start so we have infinite precision.
  size_t precision = std::numeric_limits<size_t>::max();
  Padic<T> x{eff_valuation, precision, coefficients};
  //  std::cerr << "Padic_from_integer : ";
  //  Padic_debug_print(x, std::cerr);
  return x;
}

template <typename T>
Padic<T> Padic_from_integer(T const &val, T const &p, size_t const &precision) {
  std::vector<T> coefficients;
  T val_work = val;
  int eff_valuation = 0;
  bool has_nz = false;
  while (true) {
    if (val_work == 0) {
      break;
    }
    std::pair<T, T> pair = ResQuoInt(val_work, p);
    if (pair.first == 0) {
      if (!has_nz) {
        eff_valuation += 1;
      } else {
        coefficients.push_back(pair.first);
      }
    } else {
      has_nz = true;
      coefficients.push_back(pair.first);
    }
    if (coefficients.size() == precision) {
      break;
    }
    val_work = pair.second;
  }
  return {eff_valuation, precision, coefficients};
}

template <typename T>
Padic<T> Padic_reduce_precision(Padic<T> const &x,
                                size_t const &new_precision) {
#ifdef DEBUG_PADIC
  if (x.precision < new_precision) {
    std::cerr << "The precision can only be decreased\n";
    throw TerminalException{1};
  }
#endif
  size_t len = std::min(x.coefficients.size(), new_precision);
  std::vector<T> coefficients;
  for (size_t u = 0; u < len; u++) {
    coefficients.push_back(x.coefficients[u]);
  }
  return {x.eff_valuation, new_precision, coefficients};
}

template <typename T> int Padic_valuation(Padic<T> const &x) {
  size_t infinite_precision = std::numeric_limits<size_t>::max();
  if (x.precision == infinite_precision) {
    int infinite_valuation = std::numeric_limits<int>::max();
    return infinite_valuation;
  }
  int valuation = x.eff_valuation;
  for (size_t u = 0; u < x.coefficients.size(); u++) {
    if (x.coefficients[u] == 0) {
      valuation += 1;
    } else {
      return valuation;
    }
  }
  throw PadicPrecisionException{1};
}

template <typename T> T Padic_coeff(Padic<T> const &x, size_t const &index) {
  if (index >= x.precision) {
    throw PadicPrecisionException{2};
  }
  size_t len = x.coefficients.size();
  if (index < len) {
    return x.coefficients[index];
  }
  return T(0);
}

template <typename T> Padic<T> Padic_reduction(Padic<T> const &x) {
  size_t infinite_precision = std::numeric_limits<size_t>::max();
  size_t len = x.coefficients.size();
  for (size_t u = 0; u < len; u++) {
    if (x.coefficients[u] != 0) {
      size_t precision = infinite_precision;
      if (x.precision < infinite_precision) {
        precision -= u;
      }
      std::vector<T> coefficients;
      for (size_t v = u; v < len; v++) {
        coefficients.push_back(x.coefficients[v]);
      }
      int eff_valuation = x.eff_valuation + u;
      return {eff_valuation, precision, coefficients};
    }
  }
  if (x.precision == infinite_precision) {
    // Returns an exact zero.
    return {0, infinite_precision, {}};
  } else {
    // We increase the valuation and drop the precision.
    int eff_valuation = x.eff_valuation + len;
    size_t precision = x.precision - len;
    return {eff_valuation, precision, {}};
  }
}

template <typename T>
T Padic_T_sum(Padic<T> const &x, T const &p, size_t const &precision) {
  size_t len = std::min(x.coefficients.size(), precision);
  T sum(0);
  T pow = 1;
  for (size_t u = 0; u < len; u++) {
    sum += x.coefficients[u] * pow;
    pow *= p;
  }
  return sum;
}

template <typename T>
std::vector<T> Padic_vec_T(T const &x_sum, T const &p,
                           size_t const &precision) {
  std::vector<T> coefficients(precision);
  T work_val = x_sum;
  for (size_t u = 0; u < precision; u++) {
    std::pair<T, T> pair = ResQuoInt(work_val, p);
    coefficients[u] = pair.first;
    work_val = pair.second;
  }
  return coefficients;
}

template <typename T>
Padic<T> Padic_product(Padic<T> const &x, Padic<T> const &y, T const &p) {
  size_t precision = std::min(x.precision, y.precision);
  int eff_valuation = x.eff_valuation + y.eff_valuation;
  T x_sum = Padic_T_sum(x, p, precision);
  T y_sum = Padic_T_sum(y, p, precision);
  T xy_sum = x_sum * y_sum;
  std::vector<T> coefficients = Padic_vec_T(xy_sum, p, precision);
  return {eff_valuation, precision, coefficients};
}

template <typename T> Padic<T> Padic_inverse(Padic<T> const &x, T const &p) {
  if (x.coefficients[0] == 0) {
    throw PadicPrecisionException{3};
  }
  size_t infinite_precision = std::numeric_limits<size_t>::max();
  size_t precision = x.precision;
  if (precision == infinite_precision) {
    throw PadicPrecisionException{4};
  }
  T sum = Padic_T_sum(x, p, precision);
  // full pow
  T full_pow = 1;
  for (size_t u = 0; u < x.precision; u++) {
    full_pow *= p;
  }
  T inv_val = mod_inv(sum, full_pow);
  std::vector<T> coefficients = Padic_vec_T(inv_val, p, precision);
  int eff_valuation = -x.eff_valuation;
  return {eff_valuation, precision, coefficients};
}

template <typename T>
Padic<T> Padic_addition([[maybe_unused]] Padic<T> const &x,
                        [[maybe_unused]] Padic<T> const &y,
                        [[maybe_unused]] T const &p) {
  std::cerr << "Not yet implemented\n";
  throw TerminalException{1};
}

// We assume that x is reduced
template <typename T> bool Padic_is_square(Padic<T> const &x, T const &p) {
#ifdef DEBUG_PADIC
  if (x.coefficients[0] == 0) {
    std::cerr << "First coefficient has to be non-zero\n";
    throw TerminalException{1};
  }
#endif
  int res = ResInt(x.eff_valuation, 2);
  if (res == 1) {
    return false;
  }
  //
  T two(2);
  if (p == two) {
    if (x.precision < 3) {
      throw PadicPrecisionException{5};
    }
#ifdef DEBUG_PADIC
    T coeff0 = Padic_coeff(x, 0);
    if (coeff0 != 1) {
      std::cerr << "Inconsistent value\n";
      throw TerminalException{1};
    }
#endif
    T coeff1 = Padic_coeff(x, 1);
    T coeff2 = Padic_coeff(x, 2);
    if (coeff1 != 0 || coeff2 != 0) {
      return false;
    }
    return true;
  } else {
    T coeff0 = Padic_coeff(x, 0);
    return is_quadratic_residue(coeff0, p);
  }
}

template<typename T>
std::vector<T> Padic_get_residue_classes(T const& p) {
  T two(2);
  std::vector<T> classes;
  if (p == 2) {
    std::vector<int> V{1, 3, 5, 7};
    for (auto & val_i : V) {
      T val1(val_i);
      T val2 = p * val1;
      classes.push_back(val1);
      classes.push_back(val2);
    }
  } else {
    auto get_non_residue=[&]() -> T {
      T a(2);
      while(true) {
        bool test = is_quadratic_residue(a, p);
        if (!test) {
          return a;
        }
        a += 1;
      }
    };
    std::vector<T> V{T(1), get_non_residue()};
    for (auto & val : V) {
      T val1 = val;
      T val2 = p * val;
      classes.push_back(val1);
      classes.push_back(val2);
    }
  }
  return classes;
}

// Separate the exponent from the residue.
template<typename T>
std::pair<size_t, T> Padic_decompose(T const& val, T const& p) {
  if (val == 0) {
    std::cerr << "val = 0 so we cannot decompose as p^m x\n";
    throw TerminalException{1};
  }
  T val_work = val;
  size_t expo = 0;
  while (true) {
    std::pair<T, T> pair = ResQuoInt(val_work, p);
    if (pair.first == 0) {
      expo += 1;
      val_work = pair.second;
    } else {
      return {expo, val_work};
    }
  }
}

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYPADIC_H_
// clang-format on
