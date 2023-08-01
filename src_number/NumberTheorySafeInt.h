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

// is an implementation of Z

template <> struct is_implementation_of_Z<SafeInt64> {
  static const bool value = true;
};

// is an implementation of Q

template <> struct is_implementation_of_Q<SafeInt64> {
  static const bool value = false;
};

// is_euclidean_domain property

template <> struct is_euclidean_domain<SafeInt64> {
  static const bool value = true;
};

// is_ring_field (i.e. non-zero elements are invertible)

template <> struct is_ring_field<SafeInt64> {
  static const bool value = false;
};

// is_totally_ordered (i.e. not a complex field or ring)

template <> struct is_totally_ordered<SafeInt64> {
  static const bool value = true;
};

// is exact arithmetic
// Yes, SafeInt64 is exact arithmetic since if some overflow happens
// an exception occurs and never a bad result is propagated

template <> struct is_exact_arithmetic<SafeInt64> {
  static const bool value = true;
};

//
// Underlying ring
// For some operations, we do not need divisions
// but we need some ways to convert from one setting to another
//

template <> struct underlying_ring<Rational<SafeInt64>> {
  typedef SafeInt64 ring_type;
};

//
//

template <> struct overlying_field<SafeInt64> {
  typedef Rational<SafeInt64> field_type;
};

template <> struct overlying_field<Rational<SafeInt64>> {
  typedef Rational<SafeInt64> field_type;
};

template <> struct underlying_totally_ordered_ring<Rational<SafeInt64>> {
  typedef Rational<SafeInt64> real_type;
};

template <> struct underlying_totally_ordered_ring<SafeInt64> {
  typedef SafeInt64 real_type;
};

/*
template<>
struct underlying_totally_ordered_ring<int64_t> {
  typedef int64_t real_type;
};
*/

// hash functionality

namespace std {
template <> struct hash<SafeInt64> {
  std::size_t operator()(const mpz_class &val) const {
    unsigned long int val_uli = val.get_const_val();
    return val_uli;
  }
};
template <> struct hash<Rational<SafeInt64>> {
  std::size_t operator()(const Rational<SafeInt64> &val) const {
    int64_t const& val_den = val.get_den().get_const_val();
    int64_t const& val_num = val.get_num().get_const_val();
    size_t hash1 = std::hash<int64_t>()(val_den);
    size_t hash2 = std::hash<int64_t>()(val_num);
    return hash1 + (hash2 << 6) + (hash2 >> 2);
  }
};
// clang-format off
}  // namespace std
// clang-format on

// to_string functionality

namespace std {
std::string to_string(const SafeInt64 &val) {
  std::stringstream s;
  s << val;
  std::string converted(s.str());
  return converted;
}
std::string to_string(const Rational<SafeInt64> &val) {
  int64_t const& val_den = val.get_den().get_const_val();
  int64_t const& val_num = val.get_num().get_const_val();
  std::stringstream s;
  if (val_den != 1) {
    s << val_num << "/" << val_den;
  } else {
    s << val_num;
  }
  std::string converted(s.str());
  return converted;
}
// clang-format off
}  // namespace std
// clang-format on

inline void ResInt_Kernel(Rational<SafeInt64> const &a, Rational<SafeInt64> const &b,
                          Rational<SafeInt64> &res) {
  SafeInt64 const& a_den = a.get_const_den();
  SafeInt64 const& b_den = b.get_const_den();
  
  
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

inline SafeInt64 CanonicalizationUnit(SafeInt64 const &eVal) {
  if (eVal < 0)
    return -1;
  return 1;
}

inline Rational<SafeInt64> CanonicalizationUnit(Rational<SafeInt64> const &eVal) {
  if (eVal < 0)
    return -1;
  return 1;
}

inline mpq_class T_NormGen(mpq_class const &x) { return T_abs(x); }

inline mpz_class T_NormGen(mpz_class const &x) { return T_abs(x); }

inline bool IsInteger(Rational<SafeInt64> const &x) {
  int64_t const& val = x.get_const_den().get_const_val();
  return val == 1;
}

inline Rational<SafeInt64> GetDenominator(Rational<SafeInt64> const &x) {
  SafeInt64 const& eDen = x.get_const_den();
  return Rational(eDen);
}

// We need to have nbRow as input for template reasons. But it is unused in the
// symmetric case. So, pragma statement is needed to avoid a warning being
// thrown.

inline SafeInt64 GetDenominator([[maybe_unused]] SafeInt64 const &x) {
  return 1;
}

inline SafeInt64 GetDenominator_z(Rational<SafeInt64> const &x) { return x.get_const_den(); }

inline SafeInt64 GetDenominator_z([[maybe_unused]] SafeInt64 const &x) {
  return 1;
}

inline void ScalingInteger_Kernel(stc<Rational<SafeInt64>> const &x, SafeInt64 &x_ret) {
  x_ret = x.val.get_const_den();
}

inline void ScalingInteger_Kernel([[maybe_unused]] stc<SafeInt64> const &x,
                                  SafeInt64 &x_ret) {
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

// mpq_class as input

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1, Rational<SafeInt64> &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1, double &a2) {
  int64_t const& val_den = a1.get_const_den().get_const_val();
  int64_t const& val_num = a1.get_const_num().get_const_val();
  double val_den_d = static_cast<double>(val_den);
  double val_num_d = static_cast<double>(val_num);
  return val_num_d / val_den_d;
}

void Termination_rat_safeint_not_integer(stc<Rational<SafeInt64>> const &a1) {
  if (!IsInteger(a1.val)) {
    std::string str = "a1=" + std::to_string(a1.val) + " is not an integer";
    throw ConversionException{str};
  }
}

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1, SafeInt64 &a2) {
  Termination_rat_safeint_not_integer(a1);
  a2 = a1.val.get_const_num();
}

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1, int &a2) {
  Termination_rat_safeint_not_integer(a1);
  SafeInt64 a1_z = a1.val.get_const_num();
  a2 = static_cast<int>(a1_z.get_const_val());
}

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1, uint8_t &a2) {
  Termination_rat_safeint_not_integer(a1);
  SafeInt64 a1_z = a1.val.get_const_num();
  a2 = static_cast<uint8_t>(a1_z.get_const_val());
}

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1, int8_t &a2) {
  Termination_rat_safeint_not_integer(a1);
  SafeInt64 a1_z = a1.val.get_const_num();
  a2 = static_cast<int8_t>(a1_z.get_const_val());
}

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1, uint16_t &a2) {
  Termination_rat_safeint_not_integer(a1);
  SafeInt64 a1_z = a1.val.get_const_num();
  a2 = static_cast<uint16_t>(a1_z.get_const_val());
}

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1, int16_t &a2) {
  Termination_rat_safeint_not_integer(a1);
  SafeInt64 a1_z = a1.val.get_const_num();
  a2 = static_cast<int16_t>(a1_z.get_const_val());
}

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1, uint32_t &a2) {
  Termination_rat_safeint_not_integer(a1);
  SafeInt64 a1_z = a1.val.get_const_num();
  a2 = static_cast<uint32_t>(a1_z.get_const_val());
}

#ifdef __APPLE__
// On APPLE platform, the long is not the same as the int64_t
// And we cannot match that case with a
// std::enable_if<!std::is_same_v<int64_t,long>,void>
// statement because it is not a template function.
inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1, long &a2) {
  Termination_rat_safeint_not_integer(a1);
  SafeInt64 a1_z = a1.val.get_const_num();
  a2 = static_cast<long>(a1_z.get_const_val());
}

#endif

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1, int64_t &a2) {
  Termination_rat_safeint_not_integer(a1);
  SafeInt64 a1_z = a1.val.get_const_num();
  a2 = static_cast<int64_t>(a1_z.get_const_val());
}

// uint8_t as input

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, Rational<SafeInt64> &a2) {
  SafeInt64 val(a1.val);
  a2 = Rational(val);
}

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, SafeInt64 &a2) {
  a2 = long(a1.val);
}

// uint16_t as input

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, Rational<SafeInt64> &a2) {
  SafeInt64 val(a1.val);
  a2 = Rational(val);
}

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, SafeInt64 &a2) {
  a2 = SafeInt64(a1.val);
}

// uint32_t as input

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, Rational<SafeInt64> &a2) {
  SafeInt64 val(a1.val);
  a2 = Rational(val);
}

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, SafeInt64 &a2) {
  a2 = SafeInt64(a1.val);
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

#ifdef __APPLE__
// long as input (which is not the same as int64_t on APPLE platform)

inline void TYPE_CONVERSION(stc<long> const &a1, mpq_class &a2) { a2 = convert_mpz_class_int64_t(a1.val); }

inline void TYPE_CONVERSION(stc<long> const &a1, mpz_class &a2) { a2 = convert_mpz_class_int64_t(a1.val); }

#endif

// int64_t as input

inline void TYPE_CONVERSION(stc<int64_t> const &a1, mpq_class &a2) { a2 = convert_mpz_class_int64_t(a1.val); }

inline void TYPE_CONVERSION(stc<int64_t> const &a1, mpz_class &a2) { a2 = convert_mpz_class_int64_t(a1.val); }

// int32_t as input

inline void TYPE_CONVERSION(stc<int32_t> const &a1, mpq_class &a2) { a2 = a1.val; }

inline void TYPE_CONVERSION(stc<int32_t> const &a1, mpz_class &a2) { a2 = a1.val; }

// int16_t as input

inline void TYPE_CONVERSION(stc<int16_t> const &a1, mpq_class &a2) { a2 = a1.val; }

inline void TYPE_CONVERSION(stc<int16_t> const &a1, mpz_class &a2) { a2 = a1.val; }

// int8_t as input

inline void TYPE_CONVERSION(stc<int8_t> const &a1, mpq_class &a2) { a2 = a1.val; }

inline void TYPE_CONVERSION(stc<int8_t> const &a1, mpz_class &a2) { a2 = a1.val; }

// mpz_class as input

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, SafeInt64 &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, Rational<SafeInt64> &a2) {
  a2 = Rational(a1.val);
}

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, int16_t &a2) {
  int64_t const& eVal = a1.val.get_const_val();
  a2 = static_cast<int16_t>(eVal);
}

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, int32_t &a2) {
  int64_t const& eVal = a1.val.get_const_val();
  a2 = static_cast<int32_t>(eVal);
}

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, int64_t &a2) {
  int64_t const& eVal = a1.val.get_const_val();
  a2 = static_cast<int64_t>(eVal);
}

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, double &a2) {
  a2 = a1.val.get_d();
}

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, T_uint64_t &a2) {
  int64_t const& eVal = a1.val.get_const_val();
  a2 = static_cast<T_uint64_t>(eVal);
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, SafeInt64 &a2) {
  T_uint64_t const &eVal = a1.val;
  a2 = convert_mpz_class_uint64_t(eVal);
}

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1, T_uint64_t &a2) {
  Termination_mpq_not_integer(a1);
  SafeInt64 a1_z = a1.val.get_const_num();
  a2 = static_cast<T_uint64_t>(a1_z.get_const_val());
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, Rational<SafeInt64> &a2) {
  T_uint64_t const &val1 = a1.val;
  int64_t val2 = static_cast<int64_t>(val1);
  SafeInt64 val3(val2);
  a2 = Rational(val3);
}

bool universal_square_root(SafeInt64 &ret, SafeInt64 const &val) {
  int64_t const& val_i = val.get_const_val();
  if (val_i < 0) {
    return false;
  }
  double val_d = static_cast<double>(val_i);
  double ret_d = sqrt(val_d);
  int64_t ret_i = static_cast<int64_t>(round(ret_d));
  SafeInt64 ret(ret_i);
  SafeInt64 eProd = ret * ret;
  return eProd == val;
}

bool universal_square_root(Rational<SafeInt64> &ret, Rational<SafeInt64> const &val) {
  SafeInt64 const& val_num = val.get_const_num();
  SafeInt64 const& val_den = val.get_const_den();
  SafeInt64 ret_num, ret_den;
  if (!universal_square_root(ret_num, val_num))
    return false;
  if (!universal_square_root(ret_den, val_den))
    return false;
  ret = Rational(ret_num, ret_den);
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

inline void set_to_infinity(Rational<SafeInt64> &x) {
  int64_t val1 = std::numeric_limits<int64_t>::max();
  SafeInt64 val2(val1);
  x = Rational(val2);
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
inline Rational<SafeInt64> FractionalPart(Rational<SafeInt64> const &x) {
  SafeInt64 const& eNum = x.get_const_num();
  SafeInt64 const& eDen = x.get_const_den();
  SafeInt64 r(0);
  QUO_INT(eNum, eDen, r);
  return Rational(r, eDen);
}

inline Rational<SafeInt64> Floor_safe_rat(Rational<SafeInt64> const &x) {
  Rational<SafeInt64> eFrac = FractionalPart(x);
  return x - eFrac;
}

inline Rational<SafeInt64> Ceil_safe_rat(Rational<SafeInt64> const &x) {
  Rational<SafeInt64> eFrac = FractionalPart(x);
  if (eFrac == 0)
    return x;
  return 1 + x - eFrac;
}

inline void FloorInteger(Rational<SafeInt64> const &xI, Rational<SafeInt64> &xO) {
  xO = Floor_mpq(xI);
}

inline void FloorInteger(Rational<SafeInt64> const &xI, SafeInt64 &xO) {
  Rational<SafeInt64> xO_q = Floor_mpq(xI);
  xO = xO_q.get_const_num();
}

inline void FloorInteger(Rational<SafeInt64> const &xI, int &xO) {
  Rational<SafeInt64> xO_q = Floor_mpq(xI);
  xO = static_cast<int>(xO_q.get_const_num().get_const_val());
}

inline void FloorInteger(Rational<SafeInt64> const &xI, long &xO) {
  Rational<SafeInt64> xO_q = Floor_mpq(xI);
  xO = static_cast<long>(xO_q.get_const_num().get_const_val());
}

inline void CeilInteger(Rational<SafeInt64> const &xI, Rational<SafeInt64> &xO) {
  xO = Ceil_mpq(xI);
}

inline void CeilInteger(Rational<SafeInt64> const &xI, SafeInt64 &xO) {
  Rational<SafeInt64> xO_q = Ceil_mpq(xI);
  xO = xO_q.get_const_num();
}

inline void CeilInteger(Rational<SafeInt64> const &xI, int &xO) {
  Rational<SafeInt64> xO_q = Ceil_mpq(xI);
  xO = static_cast<int>(xO_q.get_const_num().get_const_val());
}

inline void CeilInteger(Rational<SafeInt64> const &xI, long &xO) {
  Rational<SafeInt64> xO_q = Ceil_mpq(xI);
  xO = static_cast<long>(xO_q.get_const_num().get_const_val());
}

// return the nearest integer to x.
// If x is of the form y + 1/2 then it returns y.
inline Rational<SafeInt64> NearestInteger_rni(Rational<SafeInt64> const &x) {
  Rational<SafeInt64> eFrac = FractionalPart(x);
  Rational<SafeInt64> eDiff1 = eFrac;
  Rational<SafeInt64> eDiff2 = 1 - eFrac;
  Rational<SafeInt64> RetVal = x - eFrac;
  if (eDiff1 <= eDiff2) {
    return RetVal;
  } else {
    return RetVal + 1;
  }
}

inline void NearestInteger(Rational<SafeInt64> const &xI, Rational<SafeInt64> &xO) {
  xO = NearestInteger_rni(xI);
}

inline void NearestInteger(Rational<SafeInt64> const &xI, SafeInt64 &xO) {
  Rational<SafeInt64> xO_q = NearestInteger_rni(xI);
  xO = xO_q.get_const_num();
}

inline void NearestInteger(int const &xI, Rational<SafeInt64> &xO) {
  SafeInt64 val(xI);
  xO = Rational(val);
}

inline void NearestInteger(long const &xI, Rational<SafeInt64> &xO) {
  SafeInt64 val(xI);
  xO = Rational(val);
  xO = xI;
}

inline void NearestInteger(Rational<SafeInt64> const &xI, int &xO) {
  Rational<SafeInt64> xO_q = NearestInteger_rni(xI);
  xO = static_cast<int>(xO_q.get_const_num().get_const_val());
}

inline void NearestInteger(Rational<SafeInt64> const &xI, long &xO) {
  Rational<SafeInt64> xO_q = NearestInteger_rni(xI);
  xO = static_cast<long>(xO_q.get_const_num().get_const_val());
}

namespace boost::serialization {

// mpq_class

template <class Archive>
inline void load(Archive &ar, Rational<SafeInt64> &val,
                 [[maybe_unused]] const unsigned int version) {
  std::string str;
  ar &make_nvp("safe_rat", str);
  std::istringstream is(str);
  is >> val;
}

template <class Archive>
inline void save(Archive &ar, Rational<SafeInt64> const &val,
                 [[maybe_unused]] const unsigned int version) {
  std::ostringstream os;
  os << val;
  std::string str = os.str();
  ar &make_nvp("safe_rat", str);
}

template <class Archive>
inline void serialize(Archive &ar, Rational<SafeInt64> &val,
                      [[maybe_unused]] const unsigned int version) {
  split_free(ar, val, version);
}

// mpz_class

template <class Archive>
inline void load(Archive &ar, SafeInt64 &val,
                 [[maybe_unused]] const unsigned int version) {
  std::string str;
  ar &make_nvp("safe_int", str);
  std::istringstream is(str);
  is >> val;
}

template <class Archive>
inline void save(Archive &ar, SafeInt64 const &val,
                 [[maybe_unused]] const unsigned int version) {
  std::ostringstream os;
  os << val;
  std::string str = os.str();
  ar &make_nvp("safe_int", str);
}

template <class Archive>
inline void serialize(Archive &ar, SafeInt64 &val, const unsigned int version) {
  split_free(ar, val, version);
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYGMP_H_
// clang-format on
