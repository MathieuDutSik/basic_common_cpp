#ifndef INCLUDE_NUMBER_THEORY_BOOST_GMP_INT
#define INCLUDE_NUMBER_THEORY_BOOST_GMP_INT



#include <boost/multiprecision/gmp.hpp>
#include "TemplateTraits.h"
#include "ExceptionEnding.h"
#include "boost_serialization.h"
#include <iostream>

template <>
struct is_boost_mpz_int<boost::multiprecision::mpz_int> {
  static const bool value = true;
};

template <>
struct is_boost_mpq_rational<boost::multiprecision::mpq_rational> {
  static const bool value = true;
};


// hash

namespace std {
  template <>
  struct hash<boost::multiprecision::mpz_int>
  {
    std::size_t operator()(const boost::multiprecision::mpz_int& val) const
    {
      std::stringstream s;
      s << val;
      std::string converted(s.str());
      return std::hash<std::string>()(converted);
    }
  };
  template <>
  struct hash<boost::multiprecision::mpq_rational>
  {
    std::size_t operator()(const boost::multiprecision::mpq_rational& val) const
    {
      std::stringstream s;
      s << val;
      std::string converted(s.str());
      return std::hash<std::string>()(converted);
    }
  };
}

// to_string

namespace std {
  std::string to_string(const boost::multiprecision::mpz_int& e_val)
  {
    std::stringstream s;
    s << e_val;
    std::string converted(s.str());
    return converted;
  };
  std::string to_string(const boost::multiprecision::mpq_rational& e_val)
  {
    std::stringstream s;
    s << e_val;
    std::string converted(s.str());
    return converted;
  };
}

// boost serialization

namespace boost::serialization {

  // boost::multiprecision::mpq_rational

  template<class Archive>
  inline void load(Archive & ar, boost::multiprecision::mpq_rational & val, [[maybe_unused]] const unsigned int version)
  {
    std::string str;
    ar & make_nvp("mpq_rational", str);
    std::istringstream is(str);
    is >> val;
  }

  template<class Archive>
  inline void save(Archive & ar, boost::multiprecision::mpq_rational const& val, [[maybe_unused]] const unsigned int version)
  {
    std::ostringstream os;
    os << val;
    std::string str=os.str();
    ar & make_nvp("mpq", str);
  }

  template<class Archive>
  inline void serialize(Archive & ar, boost::multiprecision::mpq_rational & val, [[maybe_unused]] const unsigned int version)
  {
    split_free(ar, val, version);
  }

  // boost::multiprecision::mpz_int

  template<class Archive>
  inline void load(Archive & ar, boost::multiprecision::mpz_int & val, [[maybe_unused]] const unsigned int version)
  {
    std::string str;
    ar & make_nvp("mpz_int", str);
    std::istringstream is(str);
    is >> val;
  }

  template<class Archive>
  inline void save(Archive & ar, boost::multiprecision::mpz_int const& val, [[maybe_unused]] const unsigned int version)
  {
    std::ostringstream os;
    os << val;
    std::string str=os.str();
    ar & make_nvp("mpz", str);
  }

  template<class Archive>
  inline void serialize(Archive & ar, boost::multiprecision::mpz_int & val, const unsigned int version)
  {
    split_free(ar, val, version);
  }

}



template <>
struct is_euclidean_domain<boost::multiprecision::mpz_int> {
  static const bool value = true;
};
template <>
struct is_euclidean_domain<boost::multiprecision::mpq_rational> {
  static const bool value = true;
};


template<>
struct is_exact_arithmetic<boost::multiprecision::mpz_int> {
  static const bool value = true;
};
template<>
struct is_exact_arithmetic<boost::multiprecision::mpq_rational> {
  static const bool value = true;
};


template<>
struct is_implementation_of_Z<boost::multiprecision::mpz_int> {
  static const bool value = true;
};
template<>
struct is_implementation_of_Z<boost::multiprecision::mpq_rational> {
  static const bool value = false;
};


template <>
struct is_ring_field<boost::multiprecision::mpz_int> {
  static const bool value = false;
};
template <>
struct is_ring_field<boost::multiprecision::mpq_rational> {
  static const bool value = true;
};


template <>
struct is_totally_ordered<boost::multiprecision::mpz_int> {
  static const bool value = true;
};

template <>
struct is_totally_ordered<boost::multiprecision::mpq_rational> {
  static const bool value = true;
};


template<>
struct underlying_ring<boost::multiprecision::mpz_int> {
  typedef boost::multiprecision::mpz_int ring_type;
};
template<>
struct underlying_ring<boost::multiprecision::mpq_rational> {
  typedef boost::multiprecision::mpz_int ring_type;
};


template<>
struct overlying_field<boost::multiprecision::mpz_int> {
  typedef boost::multiprecision::mpq_rational field_type;
};
template<>
struct overlying_field<boost::multiprecision::mpq_rational> {
  typedef boost::multiprecision::mpq_rational field_type;
};


template<>
struct underlying_totally_ordered_ring<boost::multiprecision::mpz_int> {
  typedef boost::multiprecision::mpz_int real_type;
};
template<>
struct underlying_totally_ordered_ring<boost::multiprecision::mpq_rational> {
  typedef boost::multiprecision::mpq_rational real_type;
};



inline boost::multiprecision::mpz_int CanonicalizationUnit(boost::multiprecision::mpz_int const& eVal)
{
  if (eVal < 0)
    return -1;
  return 1;
}
inline boost::multiprecision::mpq_rational CanonicalizationUnit(boost::multiprecision::mpq_rational const& eVal)
{
  if (eVal < 0)
    return -1;
  return 1;
}




inline void ResInt_Kernel(boost::multiprecision::mpz_int const& a, boost::multiprecision::mpz_int const& b, boost::multiprecision::mpz_int & res)
{
  using T = boost::multiprecision::mpz_int;
  T q = a / b;
  if (a < 0 && b * q != a) {
    if (b > 0)
      q--;
    else
      q++;
  }
  res = a - q * b;
}
inline boost::multiprecision::mpz_int QuoInt(boost::multiprecision::mpz_int const& a, boost::multiprecision::mpz_int const& b)
{
  using T = boost::multiprecision::mpz_int;
  T q = a / b;
  if (a < 0 && b * q != a) {
    if (b > 0)
      q--;
    else
      q++;
  }
  return q;
}



inline std::pair<boost::multiprecision::mpq_rational,boost::multiprecision::mpq_rational> ResQuoInt_kernel(boost::multiprecision::mpq_rational const& a, boost::multiprecision::mpq_rational const& b)
{
  // a = a_n / a_d
  // b = b_n / b_d
  // a = res + q * b  with 0 <= res < |b|
  // equivalent to
  // a_n / a_d = res + q * (b_n / b_d)
  // equivalent to
  // a_n * b_d = res * a_d * b_d + (q * a_d) * b_n
  using Tf= boost::multiprecision::mpq_rational;
  using T = boost::multiprecision::mpz_int;
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
inline void ResInt_Kernel(boost::multiprecision::mpq_rational const& a, boost::multiprecision::mpq_rational const& b, boost::multiprecision::mpq_rational & res)
{
  res = ResQuoInt_kernel(a, b).first;
}
inline boost::multiprecision::mpq_rational QuoInt(boost::multiprecision::mpq_rational const& a, boost::multiprecision::mpq_rational const& b)
{
  return ResQuoInt_kernel(a, b).second;
}









inline bool IsInteger(boost::multiprecision::mpq_rational const& x)
{
  boost::multiprecision::mpz_int one = 1;
  boost::multiprecision::mpz_int eDen = denominator(x);
  return eDen == one;
}

inline boost::multiprecision::mpz_int GetDenominator([[maybe_unused]] boost::multiprecision::mpz_int const& x)
{
  return 1;
}

inline boost::multiprecision::mpq_rational GetDenominator(boost::multiprecision::mpq_rational const& x)
{
  boost::multiprecision::mpz_int eDen = denominator(x);
  boost::multiprecision::mpq_rational eDen_q = eDen;
  return eDen_q;
}





inline void TYPE_CONVERSION(boost::multiprecision::mpz_int const& a1, double & a2)
{
  a2 = a1.template convert_to<double>();
  //  std::cerr << "Missing code, write here\n";
  //  throw TerminalException{1};
}
inline void TYPE_CONVERSION(boost::multiprecision::mpz_int const& a1, int & a2)
{
  a2 = a1.template convert_to<int>();
}
inline void TYPE_CONVERSION(boost::multiprecision::mpz_int const& a1, long & a2)
{
  a2 = a1.template convert_to<long>();
}
inline void TYPE_CONVERSION(int const& a1, boost::multiprecision::mpz_int & a2)
{
  a2 = a1;
}
inline void TYPE_CONVERSION(boost::multiprecision::mpz_int const& a1, boost::multiprecision::mpz_int & a2)
{
  a2 = a1;
}


inline void TYPE_CONVERSION(boost::multiprecision::mpq_rational const& a1, double & a2)
{
  a2 = a1.template convert_to<double>();
  //  std::cerr << "Missing code, write here\n";
  //  throw TerminalException{1};
}
inline void TYPE_CONVERSION(boost::multiprecision::mpq_rational const& a1, boost::multiprecision::mpz_int & a2)
{
  if (!IsInteger(a1)) {
    std::string str = "a1=" + std::to_string(a1) + " is not an integer";
    throw ConversionException{str};
  }
  a2 = numerator(a1);
}
inline void TYPE_CONVERSION(boost::multiprecision::mpq_rational const& a1, int & a2)
{
  boost::multiprecision::mpz_int a1_z;
  TYPE_CONVERSION(a1, a1_z);
  TYPE_CONVERSION(a1_z, a2);
}
inline void TYPE_CONVERSION(boost::multiprecision::mpq_rational const& a1, long & a2)
{
  boost::multiprecision::mpz_int a1_z;
  TYPE_CONVERSION(a1, a1_z);
  TYPE_CONVERSION(a1_z, a2);
}
inline void TYPE_CONVERSION(int const& a1, boost::multiprecision::mpq_rational & a2)
{
  a2 = a1;
}
inline void TYPE_CONVERSION(long const& a1, boost::multiprecision::mpq_rational & a2)
{
  a2 = a1;
}
inline void TYPE_CONVERSION(boost::multiprecision::mpz_int const& a1, boost::multiprecision::mpq_rational & a2)
{
  a2 = a1;
}
inline void TYPE_CONVERSION(boost::multiprecision::mpq_rational const& a1, boost::multiprecision::mpq_rational & a2)
{
  a2 = a1;
}




inline boost::multiprecision::mpz_int GetDenominator_z(boost::multiprecision::mpq_rational const& x)
{
  return denominator(x);
}

inline boost::multiprecision::mpz_int GetDenominator_z(boost::multiprecision::mpz_int const& x)
{
  return 1;
}



inline boost::multiprecision::mpq_rational FractionalPart(boost::multiprecision::mpq_rational const& x)
{
  using T = boost::multiprecision::mpz_int;
  using Tf= boost::multiprecision::mpq_rational;
  T x_n = numerator(x);
  T x_d = denominator(x);
  T res;
  ResInt_Kernel(x_n, x_d, res);
  Tf res_f = res;
  Tf x_df = x_d;
  Tf ret = res_f / x_df;
  return ret;
}

inline boost::multiprecision::mpq_rational Floor_mpq_rational(boost::multiprecision::mpq_rational const& x)
{
  boost::multiprecision::mpq_rational eFrac=FractionalPart(x);
  return x-eFrac;
}

inline boost::multiprecision::mpq_rational Ceil_mpq_rational(boost::multiprecision::mpq_rational const& x)
{
  boost::multiprecision::mpq_rational eFrac=FractionalPart(x);
  if (eFrac == 0)
    return x;
  return 1 + x - eFrac;
}



inline void FloorInteger(boost::multiprecision::mpq_rational const& xI, boost::multiprecision::mpq_rational & xO)
{
  xO = Floor_mpq_rational(xI);
}
inline void FloorInteger(boost::multiprecision::mpq_rational const& xI, boost::multiprecision::mpz_int & xO)
{
  xO = numerator(Floor_mpq_rational(xI));
}
inline void FloorInteger(boost::multiprecision::mpq_rational const& xI, int & xO)
{
  boost::multiprecision::mpz_int val = numerator(Floor_mpq_rational(xI));
  xO = val.template convert_to<int>();
}
inline void FloorInteger(boost::multiprecision::mpq_rational const& xI, long & xO)
{
  boost::multiprecision::mpz_int val = numerator(Floor_mpq_rational(xI));
  xO = val.template convert_to<long>();
}

inline void CeilInteger(boost::multiprecision::mpq_rational const& xI, boost::multiprecision::mpq_rational & xO)
{
  xO = Ceil_mpq_rational(xI);
}
inline void CeilInteger(boost::multiprecision::mpq_rational const& xI, boost::multiprecision::mpz_int & xO)
{
  xO = numerator(Ceil_mpq_rational(xI));
}
inline void CeilInteger(boost::multiprecision::mpq_rational const& xI, int & xO)
{
  boost::multiprecision::mpz_int val = numerator(Ceil_mpq_rational(xI));
  xO = val.template convert_to<int>();
}
inline void CeilInteger(boost::multiprecision::mpq_rational const& xI, long & xO)
{
  boost::multiprecision::mpz_int val = numerator(Ceil_mpq_rational(xI));
  xO = val.template convert_to<long>();
}

inline boost::multiprecision::mpq_rational NearestInteger_rni(boost::multiprecision::mpq_rational const& x)
{
  boost::multiprecision::mpq_rational eFrac=FractionalPart(x);
  boost::multiprecision::mpq_rational eDiff1=eFrac;
  boost::multiprecision::mpq_rational eDiff2=1-eFrac;
  boost::multiprecision::mpq_rational RetVal=x-eFrac;
  if (eDiff1 <= eDiff2) {
    return RetVal;
  } else {
    return RetVal+1;
  }
}
inline void NearestInteger(boost::multiprecision::mpq_rational const& xI, boost::multiprecision::mpq_rational & xO)
{
  xO=NearestInteger_rni(xI);
}
inline void NearestInteger(boost::multiprecision::mpq_rational const& xI, boost::multiprecision::mpz_int & xO)
{
  boost::multiprecision::mpq_rational xO_q=NearestInteger_rni(xI);
  xO = numerator(xO_q);
}







inline void set_to_infinity(boost::multiprecision::mpq_rational & x)
{
  x = std::numeric_limits<size_t>::max();
}

inline void set_to_infinity(boost::multiprecision::mpz_int & x)
{
  x = std::numeric_limits<size_t>::max();
}

inline boost::multiprecision::mpq_rational T_NormGen(boost::multiprecision::mpq_rational const& x)
{
  if (x < 0)
    return -x;
  return x;
}

inline boost::multiprecision::mpz_int T_NormGen(boost::multiprecision::mpz_int const& x)
{
  if (x < 0)
    return -x;
  return x;
}


bool universal_square_root(boost::multiprecision::mpz_int & ret, boost::multiprecision::mpz_int const& val)
{
  using T = boost::multiprecision::mpz_int;
  ret = sqrt(val);
  T eProd = ret * ret;
  return eProd == val;
}



bool universal_square_root(boost::multiprecision::mpq_rational & ret, boost::multiprecision::mpq_rational const& val)
{
  using T = boost::multiprecision::mpz_int;
  using Tf= boost::multiprecision::mpq_rational;
  T val_n = numerator(val);
  T val_d = denominator(val);
  T ret_n, ret_d;
  if (!universal_square_root(ret_n, val_n))
    return false;
  if (!universal_square_root(ret_d, val_d))
    return false;
  ret = Tf(ret_n) / Tf(ret_d);
  return true;
}


#endif
