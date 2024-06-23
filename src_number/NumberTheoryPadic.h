// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYPADIC_H_
#define SRC_NUMBER_NUMBERTHEORYPADIC_H_
// clang-format off
#include "Temp_common.h"
#include "InputOutput.h"
#include <limits>
#include <string>
// clang-format on

/*
  Computation with p-adic numbers.

  ---

  Computing with p-adic numbers.
  We have the degree and a vector of entries of length d with all coefficients between
  0 and p-1.
  Then:
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
      Expressing x_0 y_0 = 1 + p h and expanding we can effectively find the inverse.
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
    Therefore, if we admit the result that R_2 has size 8, then we have an algorithm
    for checking
    See https://math.stackexchange.com/questions/473595/characterization-of-integers-which-has-a-2-adic-square-root
  So, if we can compute in p-adic up to some level, then we can decide whether
  entries are trivial in R_p or not.

  ---

  We do not try to git the P-adic with other data types since a constructor
  Padic x(1);
  would make sense if we were always working with the same p. In that case the "p"
  would work the same way as the d in the quadratic field Q(sqrt(d)).

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
  * valuation() that will be always greater than eff_valuation of possibly returns
     std::numeric_limits<size_t>::max() if zero.
   If there is not enough coefficients to conclude then throw a PadicPrecisionException{1};
  * coeff(i) returns the coefficient. Throw a PadicPrecisionException{2} if not accessible.

  Now about implementation:
  * We need a function accuracy_reduction that reduce the number of coefficients and the accuracy
    if non-zero. For zero entries (ones with accuracy = MAX), the accuracy is not tested.
  * We need a function is_zero that tests if a number is zero or not.
  * The computation of product returns object whose accuracy is the minimum of both accuracies
    the eff_valuation is computed as the sum of both accuracies.
  * The inverse is done in the following way:
    - Same accuracy.
    - Assumes that the first coefficient is non-zero (that is eff_valuation = valuation).
      If zero, then throw a PadicPrecisionException{3}
    - in return eff_valuation = - eff_valuation
  * The sum is the most problematic.
    - We need to take the eff_valuation and compute what can be computed.

 */

template<typename T>
sruct Padic {
  int eff_valuation;
  size_t accuracy;
  std::vector<T> coefficients;
};

template<typename T>
Padic<T> PadicFromInteger(T const& val, T const& p) {
}

template<typename T>
Padic<T> PadicProduct(Padic<T> const& x, Padic<T> const& y, T const& p) {
}

template<typename T>
Padic<T> PadicInverse(Padic<T> const& x, T const& p) {
}

template<typename T>
Padic<T> PadicReduction(Padic<T> const& x) {
}


template<typename T>
bool PadicIsSquare(Padic<T> const& x, T const& p) {
  T two(2);
  
}





// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYPADIC_H_
// clang-format on
