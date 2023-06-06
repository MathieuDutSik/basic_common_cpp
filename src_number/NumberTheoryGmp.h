// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYGMP_H_
#define SRC_NUMBER_NUMBERTHEORYGMP_H_

#include "BasicNumberTypes.h"
#include "ResidueQuotient.h"
#include "Temp_common.h"
#include "TypeConversion.h"
#include "gmpxx.h"
#include "hash_functions.h"
#include <limits>
#include <string>
#include <utility>

// We need to use is_mpq_class and is_mpz_class.
// We cannot use the std::is_same_v<T,mpq_class>  or std::is_same_v<T,mpz_class>
// because this requires mpz_class / mpq_class to be known in the scope.

template <> struct is_mpq_class<mpq_class> {
  static const bool value = true;
};

template <> struct is_mpz_class<mpz_class> {
  static const bool value = true;
};

// is an implementation of Z

template <> struct is_implementation_of_Z<mpz_class> {
  static const bool value = true;
};

template <> struct is_implementation_of_Z<mpq_class> {
  static const bool value = false;
};

// is an implementation of Q

template <> struct is_implementation_of_Q<mpz_class> {
  static const bool value = false;
};

template <> struct is_implementation_of_Q<mpq_class> {
  static const bool value = true;
};

// is_euclidean_domain property

template <> struct is_euclidean_domain<mpz_class> {
  static const bool value = true;
};

template <> struct is_euclidean_domain<mpq_class> {
  static const bool value = true;
};

// is_ring_field (i.e. non-zero elements are invertible)

template <> struct is_ring_field<mpz_class> {
  static const bool value = false;
};

template <> struct is_ring_field<mpq_class> {
  static const bool value = true;
};

// is_totally_ordered (i.e. not a complex field or ring)

template <> struct is_totally_ordered<mpz_class> {
  static const bool value = true;
};

template <> struct is_totally_ordered<mpq_class> {
  static const bool value = true;
};

// is exact arithmetic
// --- This excludes floating point types such as float/double
// --- This excludes limited types such as int/long

template <> struct is_exact_arithmetic<mpz_class> {
  static const bool value = true;
};

template <> struct is_exact_arithmetic<mpq_class> {
  static const bool value = true;
};

//
// Underlying ring
// For some operations, we do not need divisions
// but we need some ways to convert from one setting to another
//

template <> struct underlying_ring<mpq_class> {
  typedef mpz_class ring_type;
};

//
//

template <> struct overlying_field<mpz_class> {
  typedef mpq_class field_type;
};

template <> struct overlying_field<mpq_class> {
  typedef mpq_class field_type;
};

template <> struct underlying_totally_ordered_ring<mpq_class> {
  typedef mpq_class real_type;
};

template <> struct underlying_totally_ordered_ring<mpz_class> {
  typedef mpz_class real_type;
};

/*
template<>
struct underlying_totally_ordered_ring<int64_t> {
  typedef int64_t real_type;
};
*/

// hash functionality

namespace std {
template <> struct hash<mpz_class> {
  std::size_t operator()(const mpz_class &val) const {
    const int method = 2;
    if constexpr (method == 1) {
      return hash_from_stream(val);
    }
    if constexpr (method == 2) {
      unsigned long int val_uli = mpz_get_ui(val.get_mpz_t());
      return val_uli;
    }
  }
};
template <> struct hash<mpq_class> {
  std::size_t operator()(const mpq_class &val) const {
    const int method = 2;
    if constexpr (method == 1) {
      return hash_from_stream(val);
    }
    if constexpr (method == 2) {
      mpz_class val_den = val.get_den();
      mpz_class val_num = val.get_num();
      size_t hash1 = std::hash<mpz_class>()(val_den);
      size_t hash2 = std::hash<mpz_class>()(val_num);
      return hash1 + (hash2 << 6) + (hash2 >> 2);
    }
  }
};
// clang-format off
}  // namespace std
// clang-format on

// to_string functionality

namespace std {
std::string to_string(const mpz_class &e_val) {
  std::stringstream s;
  s << e_val;
  std::string converted(s.str());
  return converted;
}
std::string to_string(const mpq_class &e_val) {
  std::stringstream s;
  s << e_val;
  std::string converted(s.str());
  return converted;
}
// clang-format off
}  // namespace std
// clang-format on

// As documented in section 5.6 this is done exactly as in C int
inline void ResInt_Kernel(mpz_class const &a, mpz_class const &b,
                          mpz_class &res) {
  mpz_cdiv_r(res.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  if (b > 0 && res != 0) {
    if (b < 0)
      res -= b;
    else
      res += b;
  }
}

inline void ResInt_Kernel(mpq_class const &a, mpq_class const &b,
                          mpq_class &res) {
  mpz_class a_den = a.get_den();
  mpz_class b_den = b.get_den();
  mpz_class eGcd;
  mpz_gcd(eGcd.get_mpz_t(), a_den.get_mpz_t(), b_den.get_mpz_t());
  mpz_class eLCM = a_den * b_den / eGcd;
  mpq_class aProd = a * eLCM;
  mpq_class bProd = b * eLCM;
  mpz_class a_num = aProd.get_num();
  mpz_class b_num = bProd.get_num();
  mpz_class b_num_pos;
  if (b_num < 0) {
    b_num_pos = -b_num;
  } else {
    b_num_pos = b_num;
  }
  mpz_class res_z = a_num % b_num_pos;
  while (true) {
    if (res_z >= 0 && res_z < b_num_pos)
      break;
    if (res_z < 0)
      res_z += b_num_pos;
    if (res_z >= b_num_pos)
      res_z -= b_num_pos;
  }
  res = mpq_class(res_z) / mpq_class(eLCM);
}

inline void QUO_INT(stc<mpz_class> const &a, stc<mpz_class> const &b,
                    mpz_class &q) {
  mpz_cdiv_q(q.get_mpz_t(), a.val.get_mpz_t(), b.val.get_mpz_t());
  if (b.val > 0 && b.val * q != a.val) {
    if (b.val > 0)
      q--;
    else
      q++;
  }
}

inline void QUO_INT(stc<mpq_class> const &a, stc<mpq_class> const &b,
                    mpq_class &q) {
  mpq_class res = ResInt(a.val, b.val);
  q = (a.val - res) / b.val;
}

#include "QuoIntFcts.h"

template <typename T> std::pair<T, T> ResQuoInt(T const &a, T const &b) {
  T res = ResInt(a, b);
  T TheQ = (a - res) / b;
  return {res, TheQ};
}

inline mpz_class CanonicalizationUnit(mpz_class const &eVal) {
  if (eVal < 0)
    return -1;
  return 1;
}

inline mpq_class CanonicalizationUnit(mpq_class const &eVal) {
  if (eVal < 0)
    return -1;
  return 1;
}

inline mpq_class T_NormGen(mpq_class const &x) { return T_abs(x); }

inline mpz_class T_NormGen(mpz_class const &x) { return T_abs(x); }

inline bool IsInteger(mpq_class const &x) {
  mpz_class eDen = x.get_den();
  return eDen == 1;
}

inline mpq_class GetDenominator(mpq_class const &x) {
  mpz_class eDen = x.get_den();
  mpq_class eDen_q = eDen;
  return eDen_q;
}

// We need to have nbRow as input for template reasons. But it is unused in the
// symmetric case. So, pragma statement is needed to avoid a warning being
// thrown.

inline mpz_class GetDenominator([[maybe_unused]] mpz_class const &x) {
  return 1;
}

inline mpz_class GetDenominator_z(mpq_class const &x) { return x.get_den(); }

inline mpz_class GetDenominator_z([[maybe_unused]] mpz_class const &x) {
  return 1;
}

inline void ScalingInteger_Kernel(stc<mpq_class> const &x, mpz_class &x_ret) {
  x_ret = x.val.get_den();
}

inline void ScalingInteger_Kernel([[maybe_unused]] stc<mpz_class> const &x,
                                  mpz_class &x_ret) {
  x_ret = 1;
}

template <typename T>
inline typename std::enable_if<is_mpz_class<T>::value, T>::type
KernelGcdPair(T const &a, T const &b) {
  mpz_class eGCD;
  mpz_gcd(eGCD.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return eGCD;
}

template <typename T>
inline typename std::enable_if<is_mpz_class<T>::value, T>::type
KernelLCMpair(T const &a, T const &b) {
  mpz_class eLCM;
  mpz_lcm(eLCM.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
  return eLCM;
}

inline mpq_class GetFieldElement(mpz_class const &eVal) { return eVal; }

inline mpq_class GetFieldElement(long const &eVal) { return eVal; }

// mpq_class as input

inline void TYPE_CONVERSION(stc<mpq_class> const &a1, mpq_class &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<mpq_class> const &a1, double &a2) {
  a2 = a1.val.get_d();
}

void Termination_mpq_not_integer(stc<mpq_class> const &a1) {
  if (!IsInteger(a1.val)) {
    std::string str = "a1=" + std::to_string(a1.val) + " is not an integer";
    throw ConversionException{str};
  }
}

inline void TYPE_CONVERSION(stc<mpq_class> const &a1, mpz_class &a2) {
  Termination_mpq_not_integer(a1);
  a2 = a1.val.get_num();
}

inline void TYPE_CONVERSION(stc<mpq_class> const &a1, int &a2) {
  Termination_mpq_not_integer(a1);
  mpz_class a1_z = a1.val.get_num();
  long a1_long = a1_z.get_si();
  a2 = static_cast<int>(a1_long);
}

inline void TYPE_CONVERSION(stc<mpq_class> const &a1, uint8_t &a2) {
  Termination_mpq_not_integer(a1);
  mpz_class a1_z = a1.val.get_num();
  long a1_long = a1_z.get_si();
  a2 = uint8_t(a1_long);
}

inline void TYPE_CONVERSION(stc<mpq_class> const &a1, int8_t &a2) {
  Termination_mpq_not_integer(a1);
  mpz_class a1_z = a1.val.get_num();
  long a1_long = a1_z.get_si();
  a2 = int8_t(a1_long);
}

inline void TYPE_CONVERSION(stc<mpq_class> const &a1, uint16_t &a2) {
  Termination_mpq_not_integer(a1);
  mpz_class a1_z = a1.val.get_num();
  long a1_long = a1_z.get_si();
  a2 = uint16_t(a1_long);
}

inline void TYPE_CONVERSION(stc<mpq_class> const &a1, int16_t &a2) {
  Termination_mpq_not_integer(a1);
  mpz_class a1_z = a1.val.get_num();
  long a1_long = a1_z.get_si();
  a2 = int16_t(a1_long);
}

inline void TYPE_CONVERSION(stc<mpq_class> const &a1, uint32_t &a2) {
  Termination_mpq_not_integer(a1);
  mpz_class a1_z = a1.val.get_num();
  long a1_long = a1_z.get_si();
  a2 = uint32_t(a1_long);
}

inline void TYPE_CONVERSION(stc<mpq_class> const &a1, long &a2) {
  Termination_mpq_not_integer(a1);
  mpz_class a1_z = a1.val.get_num();
  a2 = a1_z.get_si();
}

inline void TYPE_CONVERSION(stc<mpq_class> const &a1, int64_t &a2) {
  Termination_mpq_not_integer(a1);
  mpz_class a1_z = a1.val.get_num();
  a2 = a1_z.get_si();
}

/*
inline void TYPE_CONVERSION(stc<mpq_class> const& a1, int64_t & a2)
{
  if (!IsInteger(a1.val)) {
    std::string str = "a1=" + std::to_string(a1.val) + " is not an integer";
    throw ConversionException{str};
  }
  mpz_class a1_z=a1.val.get_num();
  a2 = int64_t(a1_z.get_si());
}
*/

// uint8_t as input

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, mpq_class &a2) {
  a2 = long(a1.val);
}

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, mpz_class &a2) {
  a2 = long(a1.val);
}

// uint16_t as input

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, mpq_class &a2) {
  a2 = long(a1.val);
}

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, mpz_class &a2) {
  a2 = long(a1.val);
}

// uint32_t as input

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, mpq_class &a2) {
  a2 = long(a1.val);
}

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, mpz_class &a2) {
  a2 = long(a1.val);
}

// Conversion to mpz_class

mpz_class convert_mpz_class_uint64_t(uint64_t const& val) {
  mpz_class ret = 0;
  uint64_t shift = 256;
  int32_t shift_i = 256;
  uint64_t work_val = val;
  mpz_class prod = 1;
  while(true) {
    if (work_val == 0)
      return ret;
    uint64_t res = work_val % shift;
    uint64_t q = work_val / shift;
    int32_t res_t = static_cast<int32_t>(res);
    ret += mpz_class(res_t) * prod;
    prod *= shift_i;
    work_val = q;
  }
}

mpz_class convert_mpz_class_int64_t(int64_t const& val) {
  if (val > 0) {
    return convert_mpz_class_uint64_t(val);
  }
  return convert_mpz_class_uint64_t(-val);
}

// long long as input

inline void TYPE_CONVERSION(stc<int64_t> const &a1, mpq_class &a2) { a2 = convert_mpz_class_int64_t(a1.val); }

inline void TYPE_CONVERSION(stc<int64_t> const &a1, mpz_class &a2) { a2 = convert_mpz_class_int64_t(a1.val); }

// long as input

inline void TYPE_CONVERSION(stc<int32_t> const &a1, mpq_class &a2) { a2 = a1.val; }

inline void TYPE_CONVERSION(stc<int32_t> const &a1, mpz_class &a2) { a2 = a1.val; }

// int as input

inline void TYPE_CONVERSION(stc<int16_t> const &a1, mpq_class &a2) { a2 = a1.val; }

inline void TYPE_CONVERSION(stc<int16_t> const &a1, mpz_class &a2) { a2 = a1.val; }

// short as input

inline void TYPE_CONVERSION(stc<int8_t> const &a1, mpq_class &a2) { a2 = a1.val; }

inline void TYPE_CONVERSION(stc<int8_t> const &a1, mpz_class &a2) { a2 = a1.val; }

// mpz_class as input

inline void TYPE_CONVERSION(stc<mpz_class> const &a1, mpz_class &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<mpz_class> const &a1, mpq_class &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<mpz_class> const &a1, int16_t &a2) {
  long eVal_long = a1.val.get_si();
  a2 = static_cast<int16_t>(eVal_long);
}

inline void TYPE_CONVERSION(stc<mpz_class> const &a1, int32_t &a2) {
  long eVal_long = a1.val.get_si();
  a2 = static_cast<int32_t>(eVal_long);
}

inline void TYPE_CONVERSION(stc<mpz_class> const &a1, int64_t &a2) {
  long eVal_long = a1.val.get_si();
  a2 = static_cast<int64_t>(eVal_long);
}

inline void TYPE_CONVERSION(stc<mpz_class> const &a1, double &a2) {
  a2 = a1.val.get_d();
}

inline void TYPE_CONVERSION(stc<mpz_class> const &a1, T_uint64_t &a2) {
  long eVal_long = a1.val.get_si();
  a2 = T_uint64_t(eVal_long);
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, mpz_class &a2) {
  T_uint64_t const &eVal = a1.val;
  a2 = convert_mpz_class_uint64_t(eVal);
}

inline void TYPE_CONVERSION(stc<mpq_class> const &a1, T_uint64_t &a2) {
  Termination_mpq_not_integer(a1);
  mpz_class a1_z = a1.val.get_num();
  a2 = a1_z.get_si();
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, mpq_class &a2) {
  T_uint64_t const &eVal = a1.val;
  a2 = eVal;
}

bool universal_square_root(mpz_class &ret, mpz_class const &val) {
  mpz_sqrt(ret.get_mpz_t(), val.get_mpz_t());
  mpz_class eProd = ret * ret;
  return eProd == val;
}

bool universal_square_root(mpq_class &ret, mpq_class const &val) {
  mpz_class val_num = val.get_num();
  mpz_class val_den = val.get_den();
  mpz_class ret_num, ret_den;
  if (!universal_square_root(ret_num, val_num))
    return false;
  if (!universal_square_root(ret_den, val_den))
    return false;
  ret = mpq_class(ret_num) / mpq_class(ret_den);
  return true;
}

template <typename T> std::optional<T> UniversalSquareRoot(T const &val) {
  if (val < 0)
    return {};
  T ret;
  if (!universal_square_root(ret, val))
    return {};
  return ret;
}

inline void set_to_infinity(mpz_class &x) {
  x = std::numeric_limits<size_t>::max();
}

inline void set_to_infinity(mpq_class &x) {
  x = std::numeric_limits<size_t>::max();
}

template <typename T>
inline typename std::enable_if<std::is_integral<T>::value, void>::type
set_to_infinity(T &x) {
  x = std::numeric_limits<T>::max();
}

template <typename T> T practical_infinity() {
  T ret;
  set_to_infinity(ret);
  return ret;
}

//
// Nearest integer and similar stuff.
//
inline mpq_class FractionalPart(mpq_class const &x) {
  mpz_class eNum = x.get_num();
  mpz_class eDen = x.get_den();
  mpz_class res;
  mpz_mod(res.get_mpz_t(), eNum.get_mpz_t(), eDen.get_mpz_t());
  mpq_class eRet = mpq_class(res, eDen);
  return eRet;
}

inline mpq_class Floor_mpq(mpq_class const &x) {
  mpq_class eFrac = FractionalPart(x);
  return x - eFrac;
}

inline mpq_class Ceil_mpq(mpq_class const &x) {
  mpq_class eFrac = FractionalPart(x);
  if (eFrac == 0)
    return x;
  return 1 + x - eFrac;
}

inline void FloorInteger(mpq_class const &xI, mpq_class &xO) {
  xO = Floor_mpq(xI);
}

inline void FloorInteger(mpq_class const &xI, mpz_class &xO) {
  mpq_class xO_q = Floor_mpq(xI);
  xO = xO_q.get_num();
}

inline void FloorInteger(mpq_class const &xI, int &xO) {
  mpq_class xO_q = Floor_mpq(xI);
  xO = static_cast<int>(xO_q.get_num().get_si());
}

inline void FloorInteger(mpq_class const &xI, long &xO) {
  mpq_class xO_q = Floor_mpq(xI);
  xO = xO_q.get_num().get_si();
}

inline void CeilInteger(mpq_class const &xI, mpq_class &xO) {
  xO = Ceil_mpq(xI);
}

inline void CeilInteger(mpq_class const &xI, mpz_class &xO) {
  mpq_class xO_q = Ceil_mpq(xI);
  xO = xO_q.get_num();
}

inline void CeilInteger(mpq_class const &xI, int &xO) {
  mpq_class xO_q = Ceil_mpq(xI);
  xO = static_cast<int>(xO_q.get_num().get_si());
}

inline void CeilInteger(mpq_class const &xI, long &xO) {
  mpq_class xO_q = Ceil_mpq(xI);
  xO = xO_q.get_num().get_si();
}

// return the nearest integer to x.
// If x is of the form y + 1/2 then it returns y.
inline mpq_class NearestInteger_rni(mpq_class const &x) {
  mpq_class eFrac = FractionalPart(x);
  mpq_class eDiff1 = eFrac;
  mpq_class eDiff2 = 1 - eFrac;
  mpq_class RetVal = x - eFrac;
  if (eDiff1 <= eDiff2) {
    return RetVal;
  } else {
    return RetVal + 1;
  }
}

inline void NearestInteger(mpq_class const &xI, mpq_class &xO) {
  xO = NearestInteger_rni(xI);
}

inline void NearestInteger(mpq_class const &xI, mpz_class &xO) {
  mpq_class xO_q = NearestInteger_rni(xI);
  xO = xO_q.get_num();
}

inline void NearestInteger(int const &xI, mpq_class &xO) { xO = xI; }

inline void NearestInteger(long const &xI, mpq_class &xO) { xO = xI; }

inline void NearestInteger(mpq_class const &xI, int &xO) {
  mpq_class xO_q = NearestInteger_rni(xI);
  xO = static_cast<int>(xO_q.get_num().get_si());
}

inline void NearestInteger(mpq_class const &xI, long &xO) {
  mpq_class xO_q = NearestInteger_rni(xI);
  xO = xO_q.get_num().get_si();
}

// return the nearest integer to x.
// If x is of the form y + 1/2 then it returns y+1 and not y.
// rpi: "Rounding towards Positive Integers"
// See https://en.wikipedia.org/wiki/Floor_and_ceiling_functions#Rounding
inline mpq_class NearestInteger_rpi(mpq_class const &x) {
  //  std::cerr << "--------------------------------------\n";
  mpq_class eFrac = FractionalPart(x);
  mpq_class eOne = 1;
  mpq_class eTwo = 2;
  //  std::cerr << "We have eOne, eTwo\n";
  mpq_class eHalf = eOne / eTwo;
  //  std::cerr << "We have eHalf\n";
  mpq_class x2 = x + eHalf;
  //  std::cerr << "We have x=" << x << " eHalf=" << eHalf << " x2=" << x2 <<
  //  "\n";
  mpq_class x3 = Floor_mpq(x2);
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

template <class Archive>
inline void load(Archive &ar, mpq_class &val,
                 [[maybe_unused]] const unsigned int version) {
  std::string str;
  ar &make_nvp("mpq", str);
  std::istringstream is(str);
  is >> val;
}

template <class Archive>
inline void save(Archive &ar, mpq_class const &val,
                 [[maybe_unused]] const unsigned int version) {
  std::ostringstream os;
  os << val;
  std::string str = os.str();
  ar &make_nvp("mpq", str);
}

template <class Archive>
inline void serialize(Archive &ar, mpq_class &val,
                      [[maybe_unused]] const unsigned int version) {
  split_free(ar, val, version);
}

// mpz_class

template <class Archive>
inline void load(Archive &ar, mpz_class &val,
                 [[maybe_unused]] const unsigned int version) {
  std::string str;
  ar &make_nvp("mpz", str);
  std::istringstream is(str);
  is >> val;
}

template <class Archive>
inline void save(Archive &ar, mpz_class const &val,
                 [[maybe_unused]] const unsigned int version) {
  std::ostringstream os;
  os << val;
  std::string str = os.str();
  ar &make_nvp("mpz", str);
}

template <class Archive>
inline void serialize(Archive &ar, mpz_class &val, const unsigned int version) {
  split_free(ar, val, version);
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYGMP_H_
// clang-format on
