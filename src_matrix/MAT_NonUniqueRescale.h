// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_NONUNIQUERESCALE_H_
#define SRC_MATRIX_MAT_NONUNIQUERESCALE_H_

// clang-format off
#include "MAT_MatrixFund.h"
#include <utility>
// clang-format on

template <typename T> struct FractionVector {
  T TheMult;
  MyVector<T> TheVect;
};

template <typename T>
FractionVector<T> RemoveFractionVectorPlusCoeff(MyVector<T> const &V) {
  int n = V.size();
  auto is_zero=[&]() -> bool {
    for (int i=0; i<n; i++) {
      if (V(i) != 0) {
        return false;
      }
    }
    return true;
  };
  if (is_zero()) {
    T TheMult(1);
    return {TheMult, V};
  }
  std::vector<T> eVect(n);
  T eLCM = GetDenominator(V(0));
  for (int i = 1; i < n; i++)
    eLCM = LCMpair(eLCM, GetDenominator(V(i)));
  MyVector<T> V1 = eLCM * V;
  T eGCD = V1(0);
  for (int i = 1; i < n; i++)
    eGCD = GcdPair(eGCD, V1(i));
  MyVector<T> Vret = V1 / eGCD;
  T TheMult = eLCM / eGCD;
  return {TheMult, std::move(Vret)};
}

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
