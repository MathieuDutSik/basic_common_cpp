// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_NONUNIQUERESCALE_H_
#define SRC_MATRIX_MAT_NONUNIQUERESCALE_H_

// clang-format off
#include "MAT_MatrixFund.h"
#include <utility>
#include <vector>
// clang-format on

template <typename T> struct FractionVector {
  T TheMult;
  MyVector<T> TheVect;
};

// The scaled vector expressed over the underlying ring: the scaled entries
// are ring elements and TheVect = TheMult * V. Most uses of the scaling
// convert the scaled vector to the ring anyway, so this is the natural
// type for the result.
template <typename T> struct FractionVectorRing {
  using Tring = typename underlying_ring<T>::ring_type;
  T TheMult;
  MyVector<Tring> TheVect;
};

template <typename T>
FractionVector<T>
NonUniqueScaleToIntegerVectorPlusCoeff_Kernel(MyVector<T> const &V) {
  using Tresidual = typename T::Tresidual;
  using Tring = typename underlying_ring<Tresidual>::ring_type;
  int siz = V.size();
  Tring eLCM_ring = ScalingInteger<Tring, T>(V(0));
  for (int i = 1; i < siz; i++)
    eLCM_ring = LCMpair(eLCM_ring, ScalingInteger<Tring, T>(V(i)));
  Tresidual eLCM_res = UniversalScalarConversion<Tresidual, Tring>(eLCM_ring);
  T eLCM(eLCM_res);
  MyVector<T> Vret = V * eLCM;
  return {eLCM, std::move(Vret)};
}

// Scale the vector V to TheVect = TheMult * V with TheVect over the
// underlying ring, TheMult > 0 and, in the rational case, the entries of
// TheVect coprime (the fractions are removed and the content is reduced).
// All the arithmetic is done in the ring: no canonicalizing field
// operation occurs and the gcd of the entries is computed with an early
// exit since it usually collapses to 1 after a few entries.
template <typename T>
FractionVectorRing<T> RemoveFractionVectorPlusCoeffRing(MyVector<T> const &V) {
  using Tring = typename underlying_ring<T>::ring_type;
  // The fast branch applies to the rational types and equally to the
  // integer ring types themselves: for those there is no denominator to
  // clear and only the content reduction remains, which is what e.g.
  // CanonicalizeVector relies on.
  if constexpr (is_implementation_of_Q<T>::value ||
                is_implementation_of_Z<T>::value) {
    int n = V.size();
    if (n == 0) {
      return {T(1), MyVector<Tring>(0)};
    }
    Tring eLCM(1);
    MyVector<Tring> W(n);
    if constexpr (is_implementation_of_Z<T>::value) {
      W = UniversalVectorConversion<Tring, T>(V);
    } else {
      std::vector<Tring> dens(n);
      for (int i = 0; i < n; i++)
        dens[i] = GetDenominator_z(V(i));
      eLCM = dens[0];
      for (int i = 1; i < n; i++)
        eLCM = LCMpair(eLCM, dens[i]);
      for (int i = 0; i < n; i++)
        W(i) = GetNumerator_z(V(i)) * (eLCM / dens[i]);
    }
    Tring eGCD(0);
    for (int i = 0; i < n; i++) {
      eGCD = GcdPair(eGCD, W(i));
      if (eGCD == 1)
        break;
    }
    if constexpr (is_totally_ordered<Tring>::value) {
      if (eGCD < 0)
        eGCD = -eGCD;
    }
    if (eGCD == 0) {
      // The zero vector, for which any multiplier works.
      eGCD = 1;
    }
    if (eGCD != 1)
      for (int i = 0; i < n; i++)
        W(i) = W(i) / eGCD;
    T TheMult = UniversalScalarConversion<T, Tring>(eLCM) /
                UniversalScalarConversion<T, Tring>(eGCD);
    return {std::move(TheMult), std::move(W)};
  } else {
    FractionVector<T> fr = NonUniqueScaleToIntegerVectorPlusCoeff_Kernel(V);
    return {std::move(fr.TheMult),
            UniversalVectorConversion<Tring, T>(fr.TheVect)};
  }
}

template <typename T>
FractionVector<T> RemoveFractionVectorPlusCoeff(MyVector<T> const &V) {
  using Tring = typename underlying_ring<T>::ring_type;
  FractionVectorRing<T> frr = RemoveFractionVectorPlusCoeffRing(V);
  return {std::move(frr.TheMult),
          UniversalVectorConversion<T, Tring>(frr.TheVect)};
}

template <typename T>
FractionVector<T> NonUniqueScaleToIntegerVectorPlusCoeff(MyVector<T> const &V) {
  if constexpr (is_implementation_of_Q<T>::value) {
    return RemoveFractionVectorPlusCoeff(V);
  } else {
    return NonUniqueScaleToIntegerVectorPlusCoeff_Kernel(V);
  }
}

template <typename T>
MyVector<T> NonUniqueScaleToIntegerVector(MyVector<T> const &V) {
  return NonUniqueScaleToIntegerVectorPlusCoeff(V).TheVect;
}

// clang-format off
#endif  // SRC_MATRIX_MAT_NONUNIQUERESCALE_H_
// clang-format on
