#ifndef INCLUDE_NUMBER_GIVARO_H
#define INCLUDE_NUMBER_GIVARO_H



#include "NumberTheory.h"


#include <givaro/gfq.h>
//#include <givrational.h>
//#include <gmp++_int.h>



Givaro::Integer GetGivaroInteger(mpz_class const& val)
{
  //  mpz_t val_z(val.get_mpz_t());
  if (val >= 0) {
    size_t nchar = mpz_size(val.get_mpz_t());
    std::vector<mp_limb_t> vect_t(nchar);
    for (size_t i=0; i<nchar; i++) {
      vect_t[i] = mpz_getlimbn(val.get_mpz_t(), i);
    }
    return Givaro::Integer(vect_t);
  } else {
    mpz_class valb = -val;
    size_t nchar = mpz_size(valb.get_mpz_t());
    std::vector<mp_limb_t> vect_t(nchar);
    for (size_t i=0; i<nchar; i++) {
      vect_t[i] = mpz_getlimbn(val.get_mpz_t(), i);
    }
    return -Givaro::Integer(vect_t);
  }
}

mpz_class ConvertGivaroInteger(Givaro::Integer const& val)
{
  return mpz_class(val.get_mpz());
}





#endif
