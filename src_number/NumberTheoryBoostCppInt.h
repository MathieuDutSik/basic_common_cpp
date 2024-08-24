// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYBOOSTCPPINT_H_
#define SRC_NUMBER_NUMBERTHEORYBOOSTCPPINT_H_
// clang-format off
#include "BasicNumberTypes.h"
#include "ExceptionsFunc.h"
#include "TemplateTraits.h"
#include "TypeConversion.h"
#include "boost_serialization.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <limits>
#include <string>
#include <utility>
// clang-format on

template <> struct is_boost_cpp_int<boost::multiprecision::cpp_int> {
  static const bool value = true;
};

template <> struct is_boost_cpp_rational<boost::multiprecision::cpp_rational> {
  static const bool value = true;
};

// hash

namespace std {
template <> struct hash<boost::multiprecision::cpp_int> {
  std::size_t operator()(const boost::multiprecision::cpp_int &val) const {
    std::stringstream s;
    s << val;
    std::string converted(s.str());
    return std::hash<std::string>()(converted);
  }
};
template <> struct hash<boost::multiprecision::cpp_rational> {
  std::size_t operator()(const boost::multiprecision::cpp_rational &val) const {
    std::stringstream s;
    s << val;
    std::string converted(s.str());
    return std::hash<std::string>()(converted);
  }
};
// clang-format off
}  // namespace std
// clang-format on

// to_string

namespace std {
std::string to_string(const boost::multiprecision::cpp_int &e_val) {
  std::stringstream s;
  s << e_val;
  std::string converted(s.str());
  return converted;
}
std::string to_string(const boost::multiprecision::cpp_rational &e_val) {
  std::stringstream s;
  s << e_val;
  std::string converted(s.str());
  return converted;
}
// clang-format off
}  // namespace std
// clang-format on

// boost::serialization

/*
  That code is actually not needed.
  For some strange reasons, for cpp_rational we do not need to put the
serialization code But for mpq_rational this is apparently needed. The error
appears quite unexpectedly in the Group operator.


namespace boost::serialization {

  // boost::multiprecision::cpp_rational

  template<class Archive>
  inline void load(Archive & ar, boost::multiprecision::cpp_rational & val,
[[maybe_unused]] const unsigned int version)
  {
    std::string str;
    ar & make_nvp("cpp_rational", str);
    std::istringstream is(str);
    is >> val;
  }

  template<class Archive>
  inline void save(Archive & ar, boost::multiprecision::cpp_rational const& val,
[[maybe_unused]] const unsigned int version)
  {
    std::ostringstream os;
    os << val;
    std::string str=os.str();
    ar & make_nvp("cpp_rational", str);
  }

  template<class Archive>
  inline void serialize(Archive & ar, boost::multiprecision::cpp_rational & val,
[[maybe_unused]] const unsigned int version)
  {
    split_free(ar, val, version);
  }

  // boost::multiprecision::cpp_int

  template<class Archive>
  inline void load(Archive & ar, boost::multiprecision::cpp_int & val,
[[maybe_unused]] const unsigned int version)
  {
    std::string str;
    ar & make_nvp("cpp_int", str);
    std::istringstream is(str);
    is >> val;
  }

  template<class Archive>
  inline void save(Archive & ar, boost::multiprecision::cpp_int const& val,
[[maybe_unused]] const unsigned int version)
  {
    std::ostringstream os;
    os << val;
    std::string str=os.str();
    ar & make_nvp("mpz", str);
  }

  template<class Archive>
  inline void serialize(Archive & ar, boost::multiprecision::cpp_int & val,
const unsigned int version)
  {
    split_free(ar, val, version);
  }

}
*/

template <> struct is_euclidean_domain<boost::multiprecision::cpp_int> {
  static const bool value = true;
};
template <> struct is_euclidean_domain<boost::multiprecision::cpp_rational> {
  static const bool value = true;
};

template <> struct is_exact_arithmetic<boost::multiprecision::cpp_int> {
  static const bool value = true;
};
template <> struct is_exact_arithmetic<boost::multiprecision::cpp_rational> {
  static const bool value = true;
};

template <> struct is_implementation_of_Z<boost::multiprecision::cpp_int> {
  static const bool value = true;
};
template <> struct is_implementation_of_Z<boost::multiprecision::cpp_rational> {
  static const bool value = false;
};

template <> struct is_implementation_of_Q<boost::multiprecision::cpp_int> {
  static const bool value = false;
};
template <> struct is_implementation_of_Q<boost::multiprecision::cpp_rational> {
  static const bool value = true;
};

template <> struct is_ring_field<boost::multiprecision::cpp_int> {
  static const bool value = false;
};
template <> struct is_ring_field<boost::multiprecision::cpp_rational> {
  static const bool value = true;
};

template <> struct is_totally_ordered<boost::multiprecision::cpp_int> {
  static const bool value = true;
};

template <> struct is_totally_ordered<boost::multiprecision::cpp_rational> {
  static const bool value = true;
};

template <> struct underlying_ring<boost::multiprecision::cpp_int> {
  typedef boost::multiprecision::cpp_int ring_type;
};
template <> struct underlying_ring<boost::multiprecision::cpp_rational> {
  typedef boost::multiprecision::cpp_int ring_type;
};

template <> struct overlying_field<boost::multiprecision::cpp_int> {
  typedef boost::multiprecision::cpp_rational field_type;
};
template <> struct overlying_field<boost::multiprecision::cpp_rational> {
  typedef boost::multiprecision::cpp_rational field_type;
};

template <>
struct underlying_totally_ordered_ring<boost::multiprecision::cpp_int> {
  typedef boost::multiprecision::cpp_int real_type;
};
template <>
struct underlying_totally_ordered_ring<boost::multiprecision::cpp_rational> {
  typedef boost::multiprecision::cpp_rational real_type;
};

inline boost::multiprecision::cpp_int
CanonicalizationUnit(boost::multiprecision::cpp_int const &eVal) {
  if (eVal < 0)
    return -1;
  return 1;
}
inline boost::multiprecision::cpp_rational
CanonicalizationUnit(boost::multiprecision::cpp_rational const &eVal) {
  if (eVal < 0)
    return -1;
  return 1;
}

inline void ResInt_Kernel(boost::multiprecision::cpp_int const &a,
                          boost::multiprecision::cpp_int const &b,
                          boost::multiprecision::cpp_int &res) {
  using T = boost::multiprecision::cpp_int;
  T q = a / b;
  if (a < 0 && b * q != a) {
    if (b > 0)
      q--;
    else
      q++;
  }
  res = a - q * b;
}

void QUO_INT(stc<boost::multiprecision::cpp_int> const &a,
             stc<boost::multiprecision::cpp_int> const &b,
             boost::multiprecision::cpp_int &q) {
  q = a.val / b.val;
  if (a.val < 0 && b.val * q != a.val) {
    if (b.val > 0)
      q--;
    else
      q++;
  }
}

inline std::pair<boost::multiprecision::cpp_rational,
                 boost::multiprecision::cpp_rational>
ResQuoInt_kernel(boost::multiprecision::cpp_rational const &a,
                 boost::multiprecision::cpp_rational const &b) {
  // a = a_n / a_d
  // b = b_n / b_d
  // a = res + q * b  with 0 <= res < |b|
  // equivalent to
  // a_n / a_d = res + q * (b_n / b_d)
  // equivalent to
  // a_n * b_d = res * a_d * b_d + (q * a_d) * b_n
  using Tf = boost::multiprecision::cpp_rational;
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
  while (true) {
    if (res < 0) {
      res += b_abs;
      q -= sign;
    } else {
      if (res >= b_abs) {
        res -= b_abs;
        q += sign;
      } else {
        if (res + q * b != a) {
          std::cerr << "Error in ResQuoInt_kernel for "
                       "boost::multiprecision::cpp_rational\n";
          throw TerminalException{1};
        }
        return {res, q};
      }
    }
  }
}

inline void ResInt_Kernel(boost::multiprecision::cpp_rational const &a,
                          boost::multiprecision::cpp_rational const &b,
                          boost::multiprecision::cpp_rational &res) {
  res = ResQuoInt_kernel(a, b).first;
}

void QUO_INT(stc<boost::multiprecision::cpp_rational> const &a,
             stc<boost::multiprecision::cpp_rational> const &b,
             boost::multiprecision::cpp_rational &q) {
  q = ResQuoInt_kernel(a.val, b.val).second;
}

#include "QuoIntFcts.h"

inline bool IsInteger(boost::multiprecision::cpp_rational const &x) {
  boost::multiprecision::cpp_int one = 1;
  boost::multiprecision::cpp_int eDen = denominator(x);
  return eDen == one;
}

//

inline boost::multiprecision::cpp_int
GetNumerator([[maybe_unused]] boost::multiprecision::cpp_int const &x) {
  return x;
}

inline boost::multiprecision::cpp_rational
GetNumerator(boost::multiprecision::cpp_rational const &x) {
  boost::multiprecision::cpp_int eNum = numerator(x);
  boost::multiprecision::cpp_rational eNum_q = eNum;
  return eNum_q;
}

//

inline boost::multiprecision::cpp_int
GetDenominator([[maybe_unused]] boost::multiprecision::cpp_int const &x) {
  return 1;
}

inline boost::multiprecision::cpp_rational
GetDenominator(boost::multiprecision::cpp_rational const &x) {
  boost::multiprecision::cpp_int eDen = denominator(x);
  boost::multiprecision::cpp_rational eDen_q = eDen;
  return eDen_q;
}

//

inline boost::multiprecision::cpp_int
GetNumerator_z(boost::multiprecision::cpp_rational const &x) {
  return numerator(x);
}

inline boost::multiprecision::cpp_int
GetNumerator_z(boost::multiprecision::cpp_int const &x) {
  return x;
}

//

inline boost::multiprecision::cpp_int
GetDenominator_z(boost::multiprecision::cpp_rational const &x) {
  return denominator(x);
}

inline boost::multiprecision::cpp_int
GetDenominator_z([[maybe_unused]] boost::multiprecision::cpp_int const &x) {
  return 1;
}

//

inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_int> const &a1,
                            double &a2) {
  a2 = a1.val.template convert_to<double>();
  //  std::cerr << "Missing code, write here\n";
  //  throw TerminalException{1};
}
inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_int> const &a1,
                            int &a2) {
  a2 = a1.val.template convert_to<int>();
}
inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_int> const &a1,
                            long &a2) {
  a2 = a1.val.template convert_to<long>();
}
inline void TYPE_CONVERSION(stc<int> const &a1,
                            boost::multiprecision::cpp_int &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_int> const &a1,
                            boost::multiprecision::cpp_int &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_rational> const &a1,
                            double &a2) {
  a2 = a1.val.template convert_to<double>();
  //  std::cerr << "Missing code, write here\n";
  //  throw TerminalException{1};
}
inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_rational> const &a1,
                            boost::multiprecision::cpp_int &a2) {
  if (!IsInteger(a1.val)) {
    std::string str = "a1=" + std::to_string(a1.val) + " is not an integer";
    throw ConversionException{str};
  }
  a2 = numerator(a1.val);
}
inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_rational> const &a1,
                            int &a2) {
  boost::multiprecision::cpp_int a1_z;
  TYPE_CONVERSION(a1, a1_z);
  stc<boost::multiprecision::cpp_int> stc_a1_z{a1_z};
  TYPE_CONVERSION(stc_a1_z, a2);
}
inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_rational> const &a1,
                            long &a2) {
  boost::multiprecision::cpp_int a1_z;
  TYPE_CONVERSION(a1, a1_z);
  stc<boost::multiprecision::cpp_int> stc_a1_z{a1_z};
  TYPE_CONVERSION(stc_a1_z, a2);
}
inline void TYPE_CONVERSION(stc<int> const &a1,
                            boost::multiprecision::cpp_rational &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<long> const &a1,
                            boost::multiprecision::cpp_rational &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_int> const &a1,
                            boost::multiprecision::cpp_rational &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<boost::multiprecision::cpp_rational> const &a1,
                            boost::multiprecision::cpp_rational &a2) {
  a2 = a1.val;
}

inline void
ScalingInteger_Kernel(stc<boost::multiprecision::cpp_rational> const &x,
                      boost::multiprecision::cpp_int &x_ret) {
  x_ret = denominator(x.val);
}

inline void ScalingInteger_Kernel(
    [[maybe_unused]] stc<boost::multiprecision::cpp_int> const &x,
    boost::multiprecision::cpp_int &x_ret) {
  x_ret = 1;
}

inline boost::multiprecision::cpp_rational
FractionalPart(boost::multiprecision::cpp_rational const &x) {
  using T = boost::multiprecision::cpp_int;
  using Tf = boost::multiprecision::cpp_rational;
  T x_n = numerator(x);
  T x_d = denominator(x);
  T res;
  ResInt_Kernel(x_n, x_d, res);
  Tf res_f = res;
  Tf x_df = x_d;
  Tf ret = res_f / x_df;
  return ret;
}

inline boost::multiprecision::cpp_rational
Floor_cpp_rational(boost::multiprecision::cpp_rational const &x) {
  boost::multiprecision::cpp_rational eFrac = FractionalPart(x);
  return x - eFrac;
}

inline boost::multiprecision::cpp_rational
Ceil_cpp_rational(boost::multiprecision::cpp_rational const &x) {
  boost::multiprecision::cpp_rational eFrac = FractionalPart(x);
  if (eFrac == 0)
    return x;
  return 1 + x - eFrac;
}

inline void FloorInteger(boost::multiprecision::cpp_rational const &xI,
                         boost::multiprecision::cpp_rational &xO) {
  xO = Floor_cpp_rational(xI);
}
inline void FloorInteger(boost::multiprecision::cpp_rational const &xI,
                         boost::multiprecision::cpp_int &xO) {
  xO = numerator(Floor_cpp_rational(xI));
}
inline void FloorInteger(boost::multiprecision::cpp_rational const &xI,
                         int &xO) {
  boost::multiprecision::cpp_int val = numerator(Floor_cpp_rational(xI));
  xO = val.template convert_to<int>();
}
inline void FloorInteger(boost::multiprecision::cpp_rational const &xI,
                         long &xO) {
  boost::multiprecision::cpp_int val = numerator(Floor_cpp_rational(xI));
  xO = val.template convert_to<long>();
}

inline void CeilInteger(boost::multiprecision::cpp_rational const &xI,
                        boost::multiprecision::cpp_rational &xO) {
  xO = Ceil_cpp_rational(xI);
}
inline void CeilInteger(boost::multiprecision::cpp_rational const &xI,
                        boost::multiprecision::cpp_int &xO) {
  xO = numerator(Ceil_cpp_rational(xI));
}
inline void CeilInteger(boost::multiprecision::cpp_rational const &xI,
                        int &xO) {
  boost::multiprecision::cpp_int val = numerator(Ceil_cpp_rational(xI));
  xO = val.template convert_to<int>();
}
inline void CeilInteger(boost::multiprecision::cpp_rational const &xI,
                        long &xO) {
  boost::multiprecision::cpp_int val = numerator(Ceil_cpp_rational(xI));
  xO = val.template convert_to<long>();
}

inline boost::multiprecision::cpp_rational
NearestInteger_rni(boost::multiprecision::cpp_rational const &x) {
  boost::multiprecision::cpp_rational eFrac = FractionalPart(x);
  boost::multiprecision::cpp_rational eDiff1 = eFrac;
  boost::multiprecision::cpp_rational eDiff2 = 1 - eFrac;
  boost::multiprecision::cpp_rational RetVal = x - eFrac;
  if (eDiff1 <= eDiff2) {
    return RetVal;
  } else {
    return RetVal + 1;
  }
}
inline void NearestInteger(boost::multiprecision::cpp_rational const &xI,
                           boost::multiprecision::cpp_rational &xO) {
  xO = NearestInteger_rni(xI);
}
inline void NearestInteger(boost::multiprecision::cpp_rational const &xI,
                           boost::multiprecision::cpp_int &xO) {
  boost::multiprecision::cpp_rational xO_q = NearestInteger_rni(xI);
  xO = numerator(xO_q);
}

inline void set_to_infinity(boost::multiprecision::cpp_rational &x) {
  x = std::numeric_limits<size_t>::max();
}

inline void set_to_infinity(boost::multiprecision::cpp_int &x) {
  x = std::numeric_limits<size_t>::max();
}

inline boost::multiprecision::cpp_rational
T_NormGen(boost::multiprecision::cpp_rational const &x) {
  if (x < 0)
    return -x;
  return x;
}

inline boost::multiprecision::cpp_int
T_NormGen(boost::multiprecision::cpp_int const &x) {
  if (x < 0)
    return -x;
  return x;
}

bool universal_square_root(boost::multiprecision::cpp_int &ret,
                           boost::multiprecision::cpp_int const &val) {
  using T = boost::multiprecision::cpp_int;
  ret = sqrt(val);
  T eProd = ret * ret;
  return eProd == val;
}

bool universal_square_root(boost::multiprecision::cpp_rational &ret,
                           boost::multiprecision::cpp_rational const &val) {
  using T = boost::multiprecision::cpp_int;
  using Tf = boost::multiprecision::cpp_rational;
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

#include "NumberTheoryGeneric.h"
#include "TypeConversionFinal.h"

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYBOOSTCPPINT_H_
// clang-format on
