#ifndef INCLUDE_MPREAL_RELATED
#define INCLUDE_MPREAL_RELATED

#include "mpreal.h"
// If we are using mpreal, then somehow, we are forced to use mpq_class and mpz_class
// Those are just too related. However we included only the strictly necessary.
#include "gmpxx.h"

inline void TYPE_CONVERSION(mpq_class const& xI, mpfr::mpreal & xO)
{
  xO = mpfr::mpreal(xI.get_mpq_t());
}


void NearestInteger(mpfr::mpreal const& x, mpz_class & xO)
{
  mpfr::mpreal xRnd=round(x);
  long xRnd_ll=xRnd.toLong();
  mpz_class xRnd_z(xRnd_ll);
  auto GetErr=[&](mpz_class const& u) -> mpfr::mpreal {
    mpfr::mpreal uImg(u.get_mpz_t());
    mpfr::mpreal diff = mpfr::mpreal(uImg) - x;
    return T_abs(diff);
  };
  mpfr::mpreal err=GetErr(xRnd_z);
  while(true) {
    bool IsOK=true;
    for (int i=0; i<2; i++) {
      int shift=2*i -1;
      mpz_class xTest = xRnd_z + shift;
      mpfr::mpreal TheErr=GetErr(xRnd_z);
      if (TheErr < err) {
	IsOK=false;
	xRnd_z=xTest;
      }
    }
    if (IsOK)
      break;
  }
  xO=xRnd_z;
}



void NearestInteger(mpfr::mpreal const& x, mpq_class & x_q)
{
  mpz_class x_z;
  NearestInteger(x, x_z);
  x_q=mpq_class(x_z);
}


void NearestInteger(mpfr::mpreal const& xI, int & xO)
{
  mpz_class x_z;
  NearestInteger(xI, x_z);
  long eVal_long=x_z.get_si();
  xO=eVal_long;
}

void NearestInteger(mpfr::mpreal const& xI, long & xO)
{
  mpz_class x_z;
  NearestInteger(xI, x_z);
  xO = x_z.get_si();
}








#endif
