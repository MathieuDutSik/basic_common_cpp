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

#endif // SRC_NUMBER_ALL_FLINTXX_H_
