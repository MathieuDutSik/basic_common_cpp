// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_ALL_FLINTXX_H_
#define SRC_NUMBER_ALL_FLINTXX_H_

#include "fmpqxx.h"
#include "fmpz_vecxx.h"

/*
flint::fmpqxx operator*(int const&x, flint::fmpqxx const&y)
{
//  long int x2=x;
//  flint::fmpqxx xb=x2;
  return y*y;
}
*/

flint::fmpqxx TheProductB(int const &x, flint::fmpqxx const &y) {
  flint::fmpqxx xb = x;
  return xb * y;
}

// clang-format off
#endif  // SRC_NUMBER_ALL_FLINTXX_H_
// clang-format on
