#ifndef INCLUDE_NUMBER_THEORY
#define INCLUDE_NUMBER_THEORY

#include <iostream>
#include <sstream>

//#include "Temp_common.h"
#include "TypeConversion.h"
#include "ExceptionEnding.h"
#include <cstdlib>
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
struct overlying_field<double> {
  typedef double field_type;
};

template<>
struct overlying_field<float> {
  typedef float field_type;
};

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

/*
template<>
struct overlying_field<int64_t> {
  typedef mpq_class field_type;
};
*/

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

/*
template<>
struct underlying_totally_ordered_ring<int64_t> {
  typedef int64_t real_type;
};
*/

// hash functionality

namespace std {
  template <>
  struct hash<mpz_class>
  {
    std::size_t operator()(const mpz_class& val) const
    {
      const int method = 2;
      if constexpr(method == 1) {
        std::stringstream s;
        s << val;
        std::string converted(s.str());
        return std::hash<std::string>()(converted);
      }
      if constexpr(method == 2) {
        unsigned long int val_uli = mpz_get_ui(val.get_mpz_t());
        return val_uli;
      }
    }
  };
  template <>
  struct hash<mpq_class>
  {
    std::size_t operator()(const mpq_class& val) const
    {
      const int method = 2;
      if constexpr(method == 1) {
        std::stringstream s;
        s << val;
        std::string converted(s.str());
        return std::hash<std::string>()(converted);
      }
      if constexpr(method == 2) {
        mpz_class val_den=val.get_den();
        mpz_class val_num=val.get_num();
        size_t hash1 = std::hash<mpz_class>()(val_den);
        size_t hash2 = std::hash<mpz_class>()(val_num);
        return hash1 + (hash2 << 6) + (hash2 >> 2);
      }
    }
  };
}

// to_string functionality

namespace std {
  std::string to_string(const mpz_class& e_val)
  {
    std::stringstream s;
    s << e_val;
    std::string converted(s.str());
    return converted;
  };
  std::string to_string(const mpq_class& e_val)
  {
    std::stringstream s;
    s << e_val;
    std::string converted(s.str());
    return converted;
  };
}

// Parsing strings (It takes a std::string because we would need a const char* with a null terminated string
// and this we would not have with string_view



/*
template<>
mpz_class ParseScalar(std::string const& estr)
{
  const char* str = estr.c_str();
  mpz_class ret_val;
  int ret = mpz_set_str(ret_val.get_mpz_t(), str, 10);
  if (ret != 0) {
    std::cerr << "The ParseScalar for mpz_class failed\n";
    throw TerminalException{1};
  }
  return ret_val;
}
*/

template<typename T>
T ParseScalar(std::string const& estr)
{
  T ret_val;
  std::istringstream is(estr);
  is >> ret_val;
  return ret_val;
}

template<typename T>
void ParseScalar_inplace(std::string const& estr, T & ret_val)
{
  std::istringstream is(estr);
  is >> ret_val;
}



// The remainder and quotient of integers

template<typename T>
T ResInt_C_integer(T const& a, T const& b)
{
  T res2 = a % b;
  if (a<0 && res2 != 0) res2 += std::abs(b);
  return res2;
}




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
  T b_abs; // We cannot use std::abs which is not defined for all data types.
  if (b>0)
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

inline int ResInt(int const& a, int const& b)
{
  return ResInt_C_integer<int>(a, b);
}

inline long ResInt(long const& a, long const& b)
{
  return ResInt_C_integer<long>(a, b);
}

inline mpz_class ResInt(mpz_class const& a, mpz_class const& b)
{
  return ResInt_Generic<mpz_class>(a, b);
}


inline mpq_class ResInt(mpq_class const& a, mpq_class const& b)
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
  }
  mpq_class res_q=res_z;
  return res_q;
}


template<typename T>
T QuoInt_C_integer(T const& a, T const& b)
{
  T quo2 = a / b;
  if (a < 0 && b * quo2 != a) {
    if (b>0)
      quo2--;
    else
      quo2++;
  }
  return quo2;
}

template<typename T>
T QuoInt_Generic(T const& a, T const& b)
{
  T res = ResInt(a, b);
  return (a - res) / b;
}

inline int QuoInt(int const& a, int const& b)
{
  return QuoInt_C_integer<int>(a, b);
}

inline long QuoInt(long const& a, long const& b)
{
  return QuoInt_C_integer<long>(a, b);
}

inline mpz_class QuoInt(mpz_class const& a, mpz_class const& b)
{
  return QuoInt_Generic<mpz_class>(a, b);
}

inline mpq_class QuoInt(mpq_class const& a, mpq_class const& b)
{
  return QuoInt_Generic<mpq_class>(a, b);
}



template<typename T>
std::pair<T,T> ResQuoInt(T const& a, T const& b)
{
  T res = ResInt(a, b);
  T TheQ = (a - res) / b;
  return {res, TheQ};
}





inline int CanonicalizationUnit(int const& eVal)
{
  if (eVal < 0)
    return -1;
  return 1;
}

inline mpz_class CanonicalizationUnit(mpz_class const& eVal)
{
  if (eVal < 0)
    return -1;
  return 1;
}

inline mpq_class CanonicalizationUnit(mpq_class const& eVal)
{
  if (eVal < 0)
    return -1;
  return 1;
}



// T_Norm should always return an integer, whatever the input type
inline int T_Norm(int const& eVal)
{
  return abs(eVal);
}


inline int T_Norm(mpq_class const& x)
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

inline mpq_class T_NormGen(mpq_class const& x)
{
  return T_abs(x);
}

inline mpz_class T_NormGen(mpz_class const& x)
{
  return T_abs(x);
}

inline int T_NormGen(int const& x)
{
  return abs(x);
}







inline bool IsInteger(mpq_class const& x)
{
  mpz_class eDen=x.get_den();
  return eDen == 1;
}



inline mpq_class GetDenominator(mpq_class const& x)
{
  mpz_class eDen=x.get_den();
  mpq_class eDen_q=eDen;
  return eDen_q;
}

// We need to have nbRow as input for template reasons. But it is unused in the symmetric case.
// So, pragma statement is needed to avoid a warning being thrown.
#pragma GCC diagnostic ignored "-Wunused-parameter"
inline int GetDenominator(int const& x)
{
  return 1;
}
#pragma GCC diagnostic pop

// We need to have nbRow as input for template reasons. But it is unused in the symmetric case.
// So, pragma statement is needed to avoid a warning being thrown.
#pragma GCC diagnostic ignored "-Wunused-parameter"
inline long GetDenominator(long const& x)
{
  return 1;
}
#pragma GCC diagnostic pop

// We need to have nbRow as input for template reasons. But it is unused in the symmetric case.
// So, pragma statement is needed to avoid a warning being thrown.
#pragma GCC diagnostic ignored "-Wunused-parameter"
inline mpz_class GetDenominator(mpz_class const& x)
{
  return 1;
}
#pragma GCC diagnostic pop


inline mpz_class GetDenominator_z(mpq_class const& x)
{
  return x.get_den();
}

// We need to have nbRow as input for template reasons. But it is unused in the symmetric case.
// So, pragma statement is needed to avoid a warning being thrown.
#pragma GCC diagnostic ignored "-Wunused-parameter"
inline int GetDenominator_z(int const& x)
{
  return 1;
}
#pragma GCC diagnostic pop

// We need to have nbRow as input for template reasons. But it is unused in the symmetric case.
// So, pragma statement is needed to avoid a warning being thrown.
#pragma GCC diagnostic ignored "-Wunused-parameter"
inline long GetDenominator_z(long const& x)
{
  return 1;
}
#pragma GCC diagnostic pop

// We need to have nbRow as input for template reasons. But it is unused in the symmetric case.
// So, pragma statement is needed to avoid a warning being thrown.
#pragma GCC diagnostic ignored "-Wunused-parameter"
inline mpz_class GetDenominator_z(mpz_class const& x)
{
  return 1;
}
#pragma GCC diagnostic pop







inline mpq_class GetFieldElement(mpz_class const& eVal)
{
  return eVal;
}


inline mpq_class GetFieldElement(long const& eVal)
{
  return eVal;
}

// mpq_class as input

inline void TYPE_CONVERSION(mpq_class const& a1, mpq_class & a2)
{
  a2=a1;
}


inline void TYPE_CONVERSION(mpq_class const& a1, double & a2)
{
  a2=a1.get_d();
}


inline void TYPE_CONVERSION(mpq_class const& a1, mpz_class & a2)
{
  if (!IsInteger(a1)) {
    std::string str = "a1=" + std::to_string(a1) + " is not an integer";
    throw ConversionException{str};
  }
  a2=a1.get_num();
}

inline void TYPE_CONVERSION(mpq_class const& a1, int & a2)
{
  if (!IsInteger(a1)) {
    std::string str = "a1=" + std::to_string(a1) + " is not an integer";
    throw ConversionException{str};
  }
  mpz_class a1_z=a1.get_num();
  long a1_long=a1_z.get_si();
  a2 = int(a1_long);
}

inline void TYPE_CONVERSION(mpq_class const& a1, uint8_t & a2)
{
  if (!IsInteger(a1)) {
    std::string str = "a1=" + std::to_string(a1) + " is not an integer";
    throw ConversionException{str};
  }
  mpz_class a1_z=a1.get_num();
  long a1_long=a1_z.get_si();
  a2 = uint8_t(a1_long);
}

inline void TYPE_CONVERSION(mpq_class const& a1, int8_t & a2)
{
  if (!IsInteger(a1)) {
    std::string str = "a1=" + std::to_string(a1) + " is not an integer";
    throw ConversionException{str};
  }
  mpz_class a1_z=a1.get_num();
  long a1_long=a1_z.get_si();
  a2 = int8_t(a1_long);
}

inline void TYPE_CONVERSION(mpq_class const& a1, int16_t & a2)
{
  if (!IsInteger(a1)) {
    std::string str = "a1=" + std::to_string(a1) + " is not an integer";
    throw ConversionException{str};
  }
  mpz_class a1_z=a1.get_num();
  long a1_long=a1_z.get_si();
  a2 = int16_t(a1_long);
}

inline void TYPE_CONVERSION(mpq_class const& a1, long & a2)
{
  if (!IsInteger(a1)) {
    std::string str = "a1=" + std::to_string(a1) + " is not an integer";
    throw ConversionException{str};
  }
  mpz_class a1_z=a1.get_num();
  a2=a1_z.get_si();
}

/*
inline void TYPE_CONVERSION(mpq_class const& a1, int64_t & a2)
{
  if (!IsInteger(a1)) {
    std::string str = "a1=" + std::to_string(a1) + " is not an integer";
    throw ConversionException{str};
  }
  mpz_class a1_z=a1.get_num();
  a2 = int64_t(a1_z.get_si());
}
*/

// long as input

inline void TYPE_CONVERSION(long const& a1, long & a2)
{
  a2=a1;
}

inline void TYPE_CONVERSION(long const& a1, mpq_class & a2)
{
  a2=a1;
}

inline void TYPE_CONVERSION(long const& a1, mpz_class & a2)
{
  a2=a1;
}

// int as input

inline void TYPE_CONVERSION(int const& a1, int & a2)
{
  a2=a1;
}

inline void TYPE_CONVERSION(int const& a1, mpq_class & a2)
{
  a2=a1;
}

inline void TYPE_CONVERSION(int const& a1, mpz_class & a2)
{
  a2=a1;
}

// mpz_class as input

inline void TYPE_CONVERSION(mpz_class const& a1, mpz_class & a2)
{
  a2=a1;
}

inline void TYPE_CONVERSION(mpz_class const& a1, mpq_class & a2)
{
  a2=a1;
}

inline void TYPE_CONVERSION(mpz_class const& a1, int & a2)
{
  long eVal_long=a1.get_si();
  a2=int(eVal_long);
}

/*
inline void TYPE_CONVERSION(mpz_class const& a1, int64_t & a2)
{
  long eVal_long=a1.get_si();
  a2=int64_t(eVal_long);
}
*/

//
// Nearest integer and similar stuff.
//
inline mpq_class FractionalPart(mpq_class const& x)
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

inline mpq_class Floor(mpq_class const& x)
{
  mpq_class eFrac=FractionalPart(x);
  return x-eFrac;
}

inline mpq_class Ceil(mpq_class const& x)
{
  mpq_class eFrac=FractionalPart(x);
  if (eFrac == 0)
    return x;
  return 1 + x - eFrac;
}



// return the nearest integer to x.
// If x is of the form y + 1/2 then it returns y.
inline mpq_class NearestInteger_rni(mpq_class const& x)
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




inline void NearestInteger(mpq_class const& xI, mpq_class & xO)
{
  //  std::cerr << "NearestInteger mpq -> mpq\n";
  xO=NearestInteger_rni(xI);
}


inline void NearestInteger(mpq_class const& xI, mpz_class & xO)
{
  //  std::cerr << "NearestInteger mpq -> mpz\n";
  mpq_class xO_q=NearestInteger_rni(xI);
  xO=xO_q.get_num();
}

inline void NearestInteger(int const& xI, mpq_class & xO)
{
  xO=xI;
}

inline void NearestInteger(long const& xI, mpq_class & xO)
{
  xO=xI;
}



inline void NearestInteger(mpq_class const& xI, int & xO)
{
    mpq_class xO_q=NearestInteger_rni(xI);
    xO = int(xO_q.get_num().get_si());
}

inline void NearestInteger(mpq_class const& xI, long & xO)
{
    mpq_class xO_q=NearestInteger_rni(xI);
    xO = xO_q.get_num().get_si();
}








// return the nearest integer to x.
// If x is of the form y + 1/2 then it returns y+1 and not y.
// rpi: "Rounding towards Positive Integers"
// See https://en.wikipedia.org/wiki/Floor_and_ceiling_functions#Rounding
inline mpq_class NearestInteger_rpi(mpq_class const& x)
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





namespace boost::serialization {

  // mpq_class

  template<class Archive>
  inline void load(Archive & ar, mpq_class & val, const unsigned int version)
  {
    //      std::cerr << "load(mpq_class), step 1\n";
    std::string str;
    ar & make_nvp("mpq", str);
    std::istringstream is(str);
    is >> val;
    //      std::cerr << "load(mpq_class), step 2\n";
  }

  template<class Archive>
  inline void save(Archive & ar, mpq_class const& val, const unsigned int version)
  {
    //      std::cerr << "save(mpq_class), step 1\n";
    std::ostringstream os;
    os << val;
    std::string str=os.str();
    ar & make_nvp("mpq", str);
    //      std::cerr << "save(mpq_class), step 2\n";
  }

  template<class Archive>
  inline void serialize(Archive & ar, mpq_class & val, const unsigned int version)
  {
    //      std::cerr << "split_free(mpq_class), step 1\n";
    split_free(ar, val, version);
    //      std::cerr << "split_free(mpq_class), step 2\n";
  }

  // mpz_class

  template<class Archive>
  inline void load(Archive & ar, mpz_class & val, const unsigned int version)
  {
    //      std::cerr << "load(mpz_class), step 1\n";
    std::string str;
    ar & make_nvp("mpz", str);
    std::istringstream is(str);
    is >> val;
    //      std::cerr << "load(mpz_class), step 2\n";
  }

  template<class Archive>
  inline void save(Archive & ar, mpz_class const& val, const unsigned int version)
  {
    //      std::cerr << "save(mpz_class), step 1\n";
    std::ostringstream os;
    os << val;
    std::string str=os.str();
    ar & make_nvp("mpz", str);
    //      std::cerr << "save(mpz_class), step 2\n";
  }

  template<class Archive>
  inline void serialize(Archive & ar, mpz_class & val, const unsigned int version)
  {
    //      std::cerr << "split_free(mpz_class), step 1\n";
    split_free(ar, val, version);
    //      std::cerr << "split_free(mpz_class), step 2\n";
  }

}






#endif
