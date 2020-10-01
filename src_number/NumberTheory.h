#ifndef INCLUDE_NUMBER_THEORY
#define INCLUDE_NUMBER_THEORY

#include <iostream>

//#include "Temp_common.h"
#include "TypeConversion.h"
#include "ExceptionEnding.h"
#include "gmpxx.h"


// is an implementation of Z 

template<>
struct is_implementation_of_Z<mpz_class> {
  static const bool value = true;
};

template<>
struct is_implementation_of_Z<mpq_class> {
  static const bool value = false;
};




// is_euclidean_domain property

template <typename T>
struct is_euclidean_domain {
  static const bool value = false;
};

template <>
struct is_euclidean_domain<short> {
  static const bool value = true;
};

template <>
struct is_euclidean_domain<long> {
  static const bool value = true;
};

template <>
struct is_euclidean_domain<int> {
  static const bool value = true;
};

template <>
struct is_euclidean_domain<long long> {
  static const bool value = true;
};

template <>
struct is_euclidean_domain<mpz_class> {
  static const bool value = true;
};

template <>
struct is_euclidean_domain<mpq_class> {
  static const bool value = true;
};

// is_ring_field (i.e. non-zero elements are invertible)

template <>
struct is_ring_field<mpz_class> {
  static const bool value = false;
};

template <>
struct is_ring_field<mpq_class> {
  static const bool value = true;
};

// is_totally_ordered (i.e. not a complex field or ring)

template <>
struct is_totally_ordered<mpz_class> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<mpq_class> {
  static const bool value = true;
};

// is exact arithmetic
// --- This excludes floating point types such as float/double
// --- This excludes limited types such as int/long

template<>
struct is_exact_arithmetic<mpz_class> {
  static const bool value = true;
};

template<>
struct is_exact_arithmetic<mpq_class> {
  static const bool value = true;
};

// is_mpq_class

template <>
struct is_mpq_class<mpq_class> {
  static const bool value = true;
};

// is_mpz_class

template <>
struct is_mpz_class<mpz_class> {
  static const bool value = true;
};

//
// Underlying ring
// For some operations, we do not need divisions
// but we need some ways to convert from one setting to another
//

template<>
struct underlying_ring<mpq_class> {
  typedef mpz_class ring_type;
};


mpz_class GetRingElement(mpq_class const& eVal)
{
  return eVal.get_num();
}


//
// Overlying field
// For some operations, we do need divisions and we need
// a canonical way to do the conversion
//

template<>
struct overlying_field<mpz_class> {
  typedef mpq_class field_type;
};

template<>
struct overlying_field<int> {
  typedef mpq_class field_type;
};

template<>
struct overlying_field<long> {
  typedef mpq_class field_type;
};

template<>
struct overlying_field<mpq_class> {
  typedef mpq_class field_type;
};

template<typename T>
struct underlying_totally_ordered_ring {
};

template<>
struct underlying_totally_ordered_ring<mpq_class> {
  typedef mpq_class real_type;
};

template<>
struct underlying_totally_ordered_ring<mpz_class> {
  typedef mpz_class real_type;
};

template<>
struct underlying_totally_ordered_ring<int> {
  typedef int real_type;
};

template<>
struct underlying_totally_ordered_ring<long> {
  typedef long real_type;
};









// Specific functions for number type.
// We basically want to do the PID (Principal Ideal Domain) case.
// What are needed for that are two functions:
// ---QuotInt for the quotient
// ---T_Norm for the norm of elements.
// Result of QuoInt(a,b) is an integer q such that q is the quotient.
// We then have a = q b + r
//
// For natural integer Z (i.e. int/long/mpz_class/mpq_class)
// We should have a = bq + r
// with 0 <= r < |b| and q integer.
template<typename T>
T ResInt_Generic(T const& a, T const& b)
{
  T b_abs;
  if (b > 0)
    b_abs = b;
  else
    b_abs = -b;
  T res=a % b_abs;
  while(true) {
    if (res >= 0 && res < b_abs)
      break;
    if (res < 0)
      res += b_abs;
    if (res >= b_abs)
      res -= b_abs;
  }
  return res;
}

int ResInt(int const& a, int const& b)
{
  return ResInt_Generic<int>(a, b);
}

long ResInt(long const& a, long const& b)
{
  return ResInt_Generic<long>(a, b);
}

mpz_class ResInt(mpz_class const& a, mpz_class const& b)
{
  return ResInt_Generic<mpz_class>(a, b);
}


mpq_class ResInt(mpq_class const& a, mpq_class const& b)
{
  mpz_class a_den=a.get_den();
  mpz_class b_den=b.get_den();
  mpz_class eGcd;
  mpz_gcd(eGcd.get_mpz_t(), a_den.get_mpz_t(), b_den.get_mpz_t());
  mpz_class eLCM=a_den*b_den/eGcd;
  mpq_class aProd=a*eLCM;
  mpq_class bProd=b*eLCM;
  mpz_class a_num=aProd.get_num();
  mpz_class b_num=bProd.get_num();
  mpz_class b_num_pos;
  if (b_num < 0) {
    b_num_pos=-b_num;
  }
  else {
    b_num_pos=b_num;
  }
  mpz_class res_z=a_num % b_num_pos;
  while(true) {
    if (res_z >= 0 && res_z < b_num_pos)
      break;
    if (res_z < 0)
      res_z += b_num_pos;
    if (res_z >= b_num_pos)
      res_z -= b_num_pos;
    //    std::cerr << "res_z=" << res_z << " b_num_pos=" << b_num_pos << "\n";
  }
  //std::cerr << "a_num=" << a_num << " b_num=" << b_num << " res_z=" << res_z << "\n";
  mpq_class res_q=res_z;
  return res_q;
}

template<typename T>
T QuoInt(T const& a, T const& b)
{
  T res = ResInt(a, b);
  return (a - res) / b;
}

template<typename T>
std::pair<T,T> ResQuoInt(T const& a, T const& b)
{
  T res = ResInt(a, b);
  T TheQ = (a - res) / b;
  return {res, TheQ};
}





int CanonicalizationUnit(int const& eVal)
{
  if (eVal < 0)
    return -1;
  return 1;
}

mpz_class CanonicalizationUnit(mpz_class const& eVal)
{
  if (eVal < 0)
    return -1;
  return 1;
}

mpq_class CanonicalizationUnit(mpq_class const& eVal)
{
  if (eVal < 0)
    return -1;
  return 1;
}



// T_Norm should always return an integer, whatever the input type
int T_Norm(int const& eVal)
{
  return abs(eVal);
}


int T_Norm(mpq_class const& x)
{
  mpz_class eDen=x.get_den();
  if (eDen != 1) {
    std::cerr << "Denominator is not 1 as wished\n";
    std::cerr << "x=" << x << " eDen=" << eDen << "\n";
    std::cerr << "Error in T_Norm computation\n";
    throw TerminalException{1};
  }
  double x_d=x.get_d();
  int eValI=int(round(x_d));
  if (eValI > 0)
    return eValI;
  return -eValI;
}

mpq_class T_NormGen(mpq_class const& x)
{
  return T_abs(x);
}

mpz_class T_NormGen(mpz_class const& x)
{
  return T_abs(x);
}

int T_NormGen(int const& x)
{
  return abs(x);
}







bool IsInteger(mpq_class const& x)
{
  mpz_class eDen=x.get_den();
  return eDen == 1;
}



mpq_class GetDenominator(mpq_class const& x)
{
  mpz_class eDen=x.get_den();
  mpq_class eDen_q=eDen;
  return eDen_q;
}

int GetDenominator(int const& x)
{
  return 1;
}

long GetDenominator(long const& x)
{
  return 1;
}

mpz_class GetDenominator(mpz_class const& x)
{
  return 1;
}







mpq_class GetFieldElement(mpz_class const& eVal)
{
  return eVal;
}


mpq_class GetFieldElement(long const& eVal)
{
  return eVal;
}

// mpq_class as input

void TYPE_CONVERSION(mpq_class const& a1, mpq_class & a2)
{
  a2=a1;
}


void TYPE_CONVERSION(mpq_class const& a1, double & a2)
{
  a2=a1.get_d();
}


void TYPE_CONVERSION(mpq_class const& a1, mpz_class & a2)
{
  if (!IsInteger(a1)) {
    std::cerr << "a1=" << a1 << " is not an integer\n";
    throw TerminalException{1};
  }
  a2=a1.get_num();
}

void TYPE_CONVERSION(mpq_class const& a1, int & a2)
{
  if (!IsInteger(a1)) {
    std::cerr << "a1=" << a1 << " is not an integer\n";
    throw TerminalException{1};
  }
  mpz_class a1_z=a1.get_num();
  long a1_long=a1_z.get_si();
  a2=a1_long;
}

void TYPE_CONVERSION(mpq_class const& a1, long & a2)
{
  if (!IsInteger(a1)) {
    std::cerr << "a1=" << a1 << " is not an integer\n";
    throw TerminalException{1};
  }
  mpz_class a1_z=a1.get_num();
  a2=a1_z.get_si();
}

// long as input

void TYPE_CONVERSION(long const& a1, long & a2)
{
  a2=a1;
}

void TYPE_CONVERSION(long const& a1, mpq_class & a2)
{
  a2=a1;
}

void TYPE_CONVERSION(long const& a1, mpz_class & a2)
{
  a2=a1;
}

// int as input

void TYPE_CONVERSION(int const& a1, int & a2)
{
  a2=a1;
}

void TYPE_CONVERSION(int const& a1, mpq_class & a2)
{
  a2=a1;
}

void TYPE_CONVERSION(int const& a1, mpz_class & a2)
{
  a2=a1;
}

// mpz_class as input

void TYPE_CONVERSION(mpz_class const& a1, mpz_class & a2)
{
  a2=a1;
}

void TYPE_CONVERSION(mpz_class const& a1, mpq_class & a2)
{
  a2=a1;
}

void TYPE_CONVERSION(mpz_class const& a1, int & a2)
{
  long eVal_long=a1.get_si();
  a2=int(eVal_long);
}


//
// Nearest integer and similar stuff.
// 
mpq_class FractionalPart(mpq_class const& x)
{
  mpz_class eNum=x.get_num();
  mpz_class eDen=x.get_den();
  //  std::cerr << "FRAC eNum=" << eNum << " eDen=" << eDen << "\n";
  mpz_class res;
  mpz_mod(res.get_mpz_t(), eNum.get_mpz_t(), eDen.get_mpz_t());
  //  std::cerr << " res=" << res << "\n";
  mpq_class eRet=mpq_class(res, eDen);
  //  std::cerr << "x=" << x << " eRet=" << eRet << "\n";
  return eRet;
}

mpq_class Floor(mpq_class const& x)
{
  mpq_class eFrac=FractionalPart(x);
  return x-eFrac;
}

mpq_class Ceil(mpq_class const& x)
{
  mpq_class eFrac=FractionalPart(x);
  if (eFrac == 0)
    return x;
  return 1 + x - eFrac;
}



// return the nearest integer to x.
// If x is of the form y + 1/2 then it returns y.
mpq_class NearestInteger_rni(mpq_class const& x)
{
  mpq_class eFrac=FractionalPart(x);
  mpq_class eDiff1=eFrac;
  mpq_class eDiff2=1-eFrac;
  mpq_class RetVal=x-eFrac;
  if (eDiff1 <= eDiff2) {
    return RetVal;
  }
  else {
    return RetVal+1;
  }
}




void NearestInteger(mpq_class const& xI, mpq_class & xO)
{
  //  std::cerr << "NearestInteger mpq -> mpq\n";
  xO=NearestInteger_rni(xI);
}


void NearestInteger(mpq_class const& xI, mpz_class & xO)
{
  //  std::cerr << "NearestInteger mpq -> mpz\n";
  mpq_class xO_q=NearestInteger_rni(xI);
  xO=xO_q.get_num();
}



void NearestInteger(int const& xI, mpq_class & xO)
{
  xO=xI;
}

void NearestInteger(long const& xI, mpq_class & xO)
{
  xO=xI;
}



void NearestInteger(mpq_class const& xI, int & xO)
{
    mpq_class xO_q=NearestInteger_rni(xI);
    xO = int(xO_q.get_num().get_si());
}

void NearestInteger(mpq_class const& xI, long & xO)
{
    mpq_class xO_q=NearestInteger_rni(xI);
    xO = xO_q.get_num().get_si();
}








// return the nearest integer to x.
// If x is of the form y + 1/2 then it returns y+1 and not y.
// rpi: "Rounding towards Positive Integers"
// See https://en.wikipedia.org/wiki/Floor_and_ceiling_functions#Rounding
mpq_class NearestInteger_rpi(mpq_class const& x)
{
  //  std::cerr << "--------------------------------------\n";
  mpq_class eFrac=FractionalPart(x);
  mpq_class eOne=1;
  mpq_class eTwo=2;
  //  std::cerr << "We have eOne, eTwo\n";
  mpq_class eHalf=eOne/eTwo;
  //  std::cerr << "We have eHalf\n";
  mpq_class x2=x + eHalf;
  //  std::cerr << "We have x=" << x << " eHalf=" << eHalf << " x2=" << x2 << "\n";
  mpq_class x3=Floor(x2);
  //  std::cerr << "We have x2=" << x2 << " x3=" << x3 << "\n";
  return x3;
  /*
  mpq_class residual=x3 - eHalf;
  std::cerr << "We have x3=" << x3 << " x4=" << residual << "\n";
  std::cerr << "x=" << x << " x4=" << residual << "\n";
  if (residual < -eHalf || residual >= eHalf) {
    std::cerr << "inconsistency error in residual computation\n";
    throw TerminalException{1};
  }
  mpq_class x5=x-residual;
  std::cerr << "x=" << x << " nearestInt=" << x5 << "\n";
  std::cerr << "--------------------------------------\n";
  return x5;*/
}







#endif
