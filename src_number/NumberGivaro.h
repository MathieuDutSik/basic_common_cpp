#ifndef SRC_NUMBER_NUMBERGIVARO_H_
#define SRC_NUMBER_NUMBERGIVARO_H_

/*
  Some documentation
  https://casys.gricad-pages.univ-grenoble-alpes.fr/givaro/givaro-dev-html/class_givaro_1_1_integer.html
  https://casys.gricad-pages.univ-grenoble-alpes.fr/givaro/givaro-dev-html/class_givaro_1_1_rational.html
*/

#include "NumberTheory.h"
#include <givaro/givinteger.h>
#include <givaro/givrational.h>
#include <vector>

Givaro::Integer GetGivaroInteger(mpz_class const &val) {
  //  mpz_t val_z(val.get_mpz_t());
  auto f_convert=[&](mpz_class const& val) -> Givaro::Integer {
    size_t nchar = mpz_size(val.get_mpz_t());
    std::vector<mp_limb_t> vect_t(nchar);
    for (size_t i = 0; i < nchar; i++) {
      vect_t[i] = mpz_getlimbn(val.get_mpz_t(), i);
    }
    return Givaro::Integer(vect_t);
  };
  if (val >= 0) {
    return f_convert(val);
  } else {
    mpz_class valb = -val;
    return - f_convert(valb);
  }
}

mpz_class ConvertGivaroInteger(Givaro::Integer const &val) {
  return mpz_class(val.get_mpz());
}

Givaro::Rational GetGivaroRational(mpq_class const &val) {
  Givaro::Integer e_den = GetGivaroInteger(val.get_den());
  Givaro::Integer e_num = GetGivaroInteger(val.get_num());
  return Givaro::Rational(e_num, e_den);
}

mpq_class ConvertGivaroRational(Givaro::Rational const &val) {
  mpz_class e_num = ConvertGivaroInteger(val.nume());
  mpz_class e_den = ConvertGivaroInteger(val.deno());
  return mpq_class(e_num) / mpq_class(e_den);
}

#endif // SRC_NUMBER_NUMBERGIVARO_H_
