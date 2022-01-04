#ifndef INCLUDE_NUMBER_THEORY_BOOST
#define INCLUDE_NUMBER_THEORY_BOOST



#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/eigen.hpp>
#include "TemplateTraits.h"
#include "ExceptionEnding.h"
//#include "NumberTheory.h"



template <>
struct is_euclidean_domain<boost::multiprecision::cpp_int> {
  static const bool value = true;
};
template <>
struct is_euclidean_domain<boost::multiprecision::cpp_rational> {
  static const bool value = true;
};


template<>
struct is_exact_arithmetic<boost::multiprecision::cpp_int> {
  static const bool value = true;
};
template<>
struct is_exact_arithmetic<boost::multiprecision::cpp_rational> {
  static const bool value = true;
};


template<>
struct is_implementation_of_Z<boost::multiprecision::cpp_int> {
  static const bool value = true;
};
template<>
struct is_implementation_of_Z<boost::multiprecision::cpp_rational> {
  static const bool value = false;
};


template <>
struct is_ring_field<boost::multiprecision::cpp_int> {
  static const bool value = false;
};
template <>
struct is_ring_field<boost::multiprecision::cpp_rational> {
  static const bool value = true;
};


template <>
struct is_totally_ordered<boost::multiprecision::cpp_int> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<boost::multiprecision::cpp_rational> {
  static const bool value = true;
};


template<>
struct underlying_ring<boost::multiprecision::cpp_rational> {
  typedef boost::multiprecision::cpp_int ring_type;
};
template<>
struct overlying_field<boost::multiprecision::cpp_int> {
  typedef boost::multiprecision::cpp_rational field_type;
};


template<>
struct underlying_totally_ordered_ring<boost::multiprecision::cpp_int> {
  typedef boost::multiprecision::cpp_int real_type;
};
template<>
struct underlying_totally_ordered_ring<boost::multiprecision::cpp_rational> {
  typedef boost::multiprecision::cpp_rational real_type;
};



inline boost::multiprecision::cpp_int CanonicalizationUnit(boost::multiprecision::cpp_int const& eVal)
{
  if (eVal < 0)
    return -1;
  return 1;
}
inline boost::multiprecision::cpp_rational CanonicalizationUnit(boost::multiprecision::cpp_rational const& eVal)
{
  if (eVal < 0)
    return -1;
  return 1;
}





inline boost::multiprecision::cpp_int ResInt(boost::multiprecision::cpp_int const& a, boost::multiprecision::cpp_int const& b)
{
  using T = boost::multiprecision::cpp_int;
  T q = a / b;
  if (a < 0 && b * q != a) {
    if (b > 0)
      q--;
    else
      q++;
  }
  T res = a - q * b;
  return res;
}
inline boost::multiprecision::cpp_int QuoInt(boost::multiprecision::cpp_int const& a, boost::multiprecision::cpp_int const& b)
{
  using T = boost::multiprecision::cpp_int;
  T q = a / b;
  if (a < 0 && b * q != a) {
    if (b > 0)
      q--;
    else
      q++;
  }
  return q;
}



inline std::pair<boost::multiprecision::cpp_rational,boost::multiprecision::cpp_rational> ResQuoInt_kernel(boost::multiprecision::cpp_rational const& a, boost::multiprecision::cpp_rational const& b)
{
  // a = a_n / a_d
  // b = b_n / b_d
  // a = res + q * b  with 0 <= res < |b|
  // equivalent to
  // a_n / a_d = res + q * (b_n / b_d)
  // equivalent to
  // a_n * b_d = res * a_d * b_d + (q * a_d) * b_n
  using Tf= boost::multiprecision::cpp_rational;
  using T = boost::multiprecision::cpp_int;
  T a_n = numerator(a);
  T b_n = numerator(b);
  T a_d = denominator(a);
  T b_d = denominator(b);
  T a1 = a_n * b_d;
  T b1 = a_d * b_n;
  T q = a1 / b1;
  Tf q_f = q;
  Tf res = a - q_f * b;
  int sign;
  Tf b_abs;
  if (b < 0) {
    sign = -1;
    b_abs = -b;
  } else {
    sign = 1;
    b_abs = b;
  }
  while(true) {
    if (res < 0) {
      res += b_abs;
      q -= sign;
    } else {
      if (res >= b_abs) {
        res -= b_abs;
        q += sign;
      } else {
        if (res + q * b != a) {
          std::cerr << "Some error somewhere\n";
          throw TerminalException{1};
        }
        return {res,q};
      }
    }
  }
}
inline boost::multiprecision::cpp_rational ResInt(boost::multiprecision::cpp_rational const& a, boost::multiprecision::cpp_rational const& b)
{
  return ResQuoInt_kernel(a, b).first;
}
inline boost::multiprecision::cpp_rational QuoInt(boost::multiprecision::cpp_rational const& a, boost::multiprecision::cpp_rational const& b)
{
  return ResQuoInt_kernel(a, b).second;
}









inline bool IsInteger(boost::multiprecision::cpp_rational const& x)
{
  boost::multiprecision::cpp_int one = 1;
  boost::multiprecision::cpp_int eDen = denominator(x);
  return eDen == one;
}
inline boost::multiprecision::cpp_rational GetDenominator(boost::multiprecision::cpp_rational const& x)
{
  boost::multiprecision::cpp_int eDen = denominator(x);
  boost::multiprecision::cpp_rational eDen_q = eDen;
  return eDen_q;
}





inline void TYPE_CONVERSION(boost::multiprecision::cpp_int const& a1, double & a2)
{
  a2 = a1.template convert_to<double>();
  //  std::cerr << "Missing code, write here\n";
  //  throw TerminalException{1};
}
inline void TYPE_CONVERSION(boost::multiprecision::cpp_int const& a1, int & a2)
{
  a2 = a1.template convert_to<int>();
}
inline void TYPE_CONVERSION(boost::multiprecision::cpp_int const& a1, long & a2)
{
  a2 = a1.template convert_to<long>();
}
inline void TYPE_CONVERSION(int const& a1, boost::multiprecision::cpp_int & a2)
{
  a2 = a1;
}
inline void TYPE_CONVERSION(boost::multiprecision::cpp_int const& a1, boost::multiprecision::cpp_int & a2)
{
  a2 = a1;
}


inline void TYPE_CONVERSION(boost::multiprecision::cpp_rational const& a1, double & a2)
{
  a2 = a1.template convert_to<double>();
  //  std::cerr << "Missing code, write here\n";
  //  throw TerminalException{1};
}
inline void TYPE_CONVERSION(boost::multiprecision::cpp_rational const& a1, boost::multiprecision::cpp_int & a2)
{
  if (!IsInteger(a1)) {
    std::cerr << "a1 should be integral\n";
    throw TerminalException{1};
  }
  a2 = numerator(a1);
}
inline void TYPE_CONVERSION(boost::multiprecision::cpp_rational const& a1, int & a2)
{
  boost::multiprecision::cpp_int a1_z;
  TYPE_CONVERSION(a1, a1_z);
  TYPE_CONVERSION(a1_z, a2);
}
inline void TYPE_CONVERSION(boost::multiprecision::cpp_rational const& a1, long & a2)
{
  boost::multiprecision::cpp_int a1_z;
  TYPE_CONVERSION(a1, a1_z);
  TYPE_CONVERSION(a1_z, a2);
}
inline void TYPE_CONVERSION(int const& a1, boost::multiprecision::cpp_rational & a2)
{
  a2 = a1;
}
inline void TYPE_CONVERSION(long const& a1, boost::multiprecision::cpp_rational & a2)
{
  a2 = a1;
}
inline void TYPE_CONVERSION(boost::multiprecision::cpp_int const& a1, boost::multiprecision::cpp_rational & a2)
{
  a2 = a1;
}
inline void TYPE_CONVERSION(boost::multiprecision::cpp_rational const& a1, boost::multiprecision::cpp_rational & a2)
{
  a2 = a1;
}




inline boost::multiprecision::cpp_int GetDenominator_z(boost::multiprecision::cpp_rational const& x)
{
  return denominator(x);
}

inline boost::multiprecision::cpp_int GetDenominator_z(boost::multiprecision::cpp_int const& x)
{
  return 1;
}





#endif
