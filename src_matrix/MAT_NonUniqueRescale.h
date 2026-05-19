// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_MATRIX_MAT_NONUNIQUERESCALE_H_
#define SRC_MATRIX_MAT_NONUNIQUERESCALE_H_

// clang-format off
#include "MAT_Matrix.h"
#include <utility>
// clang-format on

template <typename T> struct FractionVector;

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
