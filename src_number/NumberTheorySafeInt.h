// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYSAFEINT_H_
#define SRC_NUMBER_NUMBERTHEORYSAFEINT_H_

// clang-format off
#include "BasicNumberTypes.h"
#include "ResidueQuotient.h"
#include "Temp_common.h"
#include "TypeConversion.h"
#include "gmpxx.h"
#include "hash_functions.h"
#include "rational.h"
#include <limits>
#include <string>
#include <utility>
// clang-format on

#define MAX_INT64_PROD 2147483647
#define MAX_INT64_SUM 4611686018427387903

void check_prod_int64(int64_t val) {
  if (val > MAX_INT64_PROD || val < -MAX_INT64_PROD) {
    std::cerr << "Safety check triggerred for product operation of int64_t\n";
    std::cerr << "val=" << val << " MAX_INT64_PROD=" << MAX_INT64_PROD << "\n";
    throw SafeIntException{val};
  }
}

void check_sum_int64(int64_t val) {
  if (val > MAX_INT64_SUM || val < -MAX_INT64_SUM) {
    std::cerr << "Safety check triggerred for sum operation of int64_t\n";
    std::cerr << "val=" << val << " MAX_INT64_SUM=" << MAX_INT64_SUM << "\n";
    throw SafeIntException{val};
  }
}

void check_reasonableness(std::string oper, int64_t val) {
  int64_t limit = MAX_INT64_SUM;
  if (T_abs(val) > limit) {
    std::cerr << "We have val=" << val << " int(val)=" << int(val) << "\n";
    std::cerr << "which does not seem reasonable at oper=" << oper << "\n";
    throw TerminalException{1};
  }
}

struct SafeInt64 {
public:
  using Tint = int64_t;

private:
  Tint val;

public:
  explicit SafeInt64() : val(0) {
#ifdef CHECK_SAFETY_INTEGER
    std::cerr << "Default constructor for SafeInt64 val=" << val << "\n";
#endif
  }
  SafeInt64(int const &x) : val(x) {
#ifdef CHECK_SAFETY_INTEGER
    std::cerr << "Int64_t constructor for SafeInt64 val=" << x << "\n";
#endif
  }
  explicit SafeInt64(Tint const &x) : val(x) {
#ifdef CHECK_SAFETY_INTEGER
    std::cerr << "Int64_t constructor for SafeInt64 val=" << x << "\n";
#endif
  }
  SafeInt64(SafeInt64 const &x) : val(x.val) {
#ifdef CHECK_SAFETY_INTEGER
    std::cerr << "SafeInt64 constructor for SafeInt64 val=" << x << "\n";
#endif
  }
  double get_d() const { return static_cast<double>(val); }
  // Assignment operators
  SafeInt64 operator=(SafeInt64 const &u) {
    // assignment operator from int
    val = u.val;
    return SafeInt64(val);
  }
  Tint &get_val() { return val; }
  const Tint &get_const_val() const { return val; }
  void operator+=(SafeInt64 const &x) {
    check_sum_int64(val);
    check_sum_int64(x.val);
    val += x.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator+", val);
#endif
  }
  void operator-=(SafeInt64 const &x) {
    check_sum_int64(val);
    check_sum_int64(x.val);
    val -= x.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator-=", val);
#endif
  }
  friend SafeInt64 operator+(SafeInt64 const &x, SafeInt64 const &y) {
    check_sum_int64(x.val);
    check_sum_int64(y.val);
    SafeInt64 z;
    z.val = x.val + y.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator+(SafeInt64,SafeInt64)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator+(Tint const &x, SafeInt64 const &y) {
    check_sum_int64(x);
    check_sum_int64(y.val);
    SafeInt64 z;
    z.val = x + y.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator+(int64_t,SafeInt64)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator+(SafeInt64 const &x, Tint const &y) {
    check_sum_int64(x.val);
    check_sum_int64(y);
    SafeInt64 z;
    z.val = x.val + y;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator+(SafeInt64,int64_t)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator-(SafeInt64 const &x, SafeInt64 const &y) {
    check_sum_int64(x.val);
    check_sum_int64(y.val);
    SafeInt64 z;
    z.val = x.val - y.val;
#ifdef CHECK_SAFETY_INTEGER
    std::cerr << "x.val=" << x.val << " y.val=" << y.val << "\n";
    check_reasonableness("operator-(SafeInt64,SafeInt64)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator-(SafeInt64 const &x) {
    SafeInt64 z;
    z.val = -x.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator-(SafeInt64)", z.val);
#endif
    return z;
  }
  SafeInt64 operator++() {
    check_sum_int64(val);
    val++;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator++()", val);
#endif
    return *this;
  }
  SafeInt64 operator++(int) {
    SafeInt64 tmp = *this;
    check_sum_int64(val);
    val++;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator++(int)", val);
#endif
    return tmp;
  }
  SafeInt64 operator--() {
    check_sum_int64(val);
    val--;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator--()", val);
#endif
    return *this;
  }
  SafeInt64 operator--(int) {
    SafeInt64 tmp = *this;
    check_sum_int64(val);
    val--;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator--(int)", val);
#endif
    return tmp;
  }
  void operator*=(SafeInt64 const &x) {
    check_sum_int64(x.val);
#ifdef CHECK_SAFETY_INTEGER
    std::cerr << "val=" << val << "\n";
#endif
    val *= x.val;
#ifdef CHECK_SAFETY_INTEGER
    std::cerr << "x.val=" << x.val << "\n";
    check_reasonableness("operator*=()", val);
#endif
  }
  friend SafeInt64 operator*(SafeInt64 const &x, SafeInt64 const &y) {
    check_prod_int64(x.val);
    check_prod_int64(y.val);
    SafeInt64 z;
    z.val = x.val * y.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator*(SafeInt64_t,SafeInt64_t)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator/(SafeInt64 const &x, SafeInt64 const &y) {
    SafeInt64 z;
    z.val = x.val / y.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator/(SafeInt64_t,SafeInt64_t)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator%(SafeInt64 const &x, SafeInt64 const &y) {
    SafeInt64 z;
    z.val = x.val % y.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator/(SafeInt64_t,SafeInt64_t)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator*(Tint const &x, SafeInt64 const &y) {
    check_prod_int64(x);
    check_prod_int64(y.val);
    SafeInt64 z;
    z.val = x * y.val;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator*(int64_t,SafeInt64)", z.val);
#endif
    return z;
  }
  friend SafeInt64 operator*(SafeInt64 const &x, Tint const &y) {
    check_prod_int64(x.val);
    check_prod_int64(y);
    SafeInt64 z;
    z.val = x.val * y;
#ifdef CHECK_SAFETY_INTEGER
    check_reasonableness("operator*(SafeInt64,int64_t)", z.val);
#endif
    return z;
  }
  friend std::ostream &operator<<(std::ostream &os, SafeInt64 const &v) {
    return os << v.val;
  }
  friend std::istream &operator>>(std::istream &is, SafeInt64 &v) {
    is >> v.val;
    return is;
  }
  //
  friend bool operator==(SafeInt64 const &x, SafeInt64 const &y) {
    return x.val == y.val;
  }
  friend bool operator==(SafeInt64 const &x, int const &y) {
    Tint y_i = static_cast<Tint>(y);
    return x.val == y_i;
  }
  //
  friend bool operator!=(SafeInt64 const &x, SafeInt64 const &y) {
    return x.val != y.val;
  }
  friend bool operator!=(SafeInt64 const &x, int const &y) {
    Tint y_i = static_cast<Tint>(y);
    return x.val != y_i;
  }
  //
  friend bool operator>=(SafeInt64 const &x, SafeInt64 const &y) {
    return x.val >= y.val;
  }
  friend bool operator>=(SafeInt64 const &x, int const &y) {
    Tint y_i = static_cast<Tint>(y);
    return x.val >= y_i;
  }
  //
  friend bool operator<=(SafeInt64 const &x, SafeInt64 const &y) {
    return x.val <= y.val;
  }
  friend bool operator<=(SafeInt64 const &x, int const &y) {
    Tint y_i = static_cast<Tint>(y);
    return x.val <= y_i;
  }
  //
  friend bool operator>(SafeInt64 const &x, SafeInt64 const &y) {
    return x.val > y.val;
  }
  friend bool operator>(SafeInt64 const &x, int const &y) {
    Tint y_i = static_cast<Tint>(y);
    return x.val > y_i;
  }
  //
  friend bool operator<(SafeInt64 const &x, SafeInt64 const &y) {
    return x.val < y.val;
  }
  friend bool operator<(SafeInt64 const &x, int const &y) {
    Tint y_i = static_cast<Tint>(y);
    return x.val < y_i;
  }
};

inline void QUO_INT(stc<SafeInt64> const &a, stc<SafeInt64> const &b,
                    SafeInt64 &q) {
  using Tint = typename SafeInt64::Tint;
  const Tint &a_int = a.val.get_const_val();
  const Tint &b_int = b.val.get_const_val();
  Tint quot = QuoInt_C_integer<Tint>(a_int, b_int);
  q = SafeInt64(quot);
}

size_t get_bit(SafeInt64 const &x) {
  int64_t const &val = x.get_const_val();
  return get_bit(val);
}

#include "QuoIntFcts.h"

void ResInt_Kernel(SafeInt64 const &a, SafeInt64 const &b, SafeInt64 &res) {
  using Tint = typename SafeInt64::Tint;
  const Tint &a_int = a.get_const_val();
  const Tint &b_int = b.get_const_val();
  res.get_val() = ResInt_C_integer<Tint>(a_int, b_int);
}

template <> struct has_reduction_subset_solver<Rational<SafeInt64>> {
  static const bool value = true;
};

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
// hash functionality
template <> struct hash<SafeInt64> {
  std::size_t operator()(const SafeInt64 &val) const {
    unsigned long int val_uli = val.get_const_val();
    return val_uli;
  }
};
template <> struct hash<Rational<SafeInt64>> {
  std::size_t operator()(const Rational<SafeInt64> &val) const {
    int64_t const &val_den = val.get_const_den().get_const_val();
    int64_t const &val_num = val.get_const_num().get_const_val();
    size_t hash1 = std::hash<int64_t>()(val_den);
    size_t hash2 = std::hash<int64_t>()(val_num);
    return hash1 + (hash2 << 6) + (hash2 >> 2);
  }
};
// to_string functionality
std::string to_string(const SafeInt64 &val) {
  std::stringstream s;
  s << val;
  std::string converted(s.str());
  return converted;
}
std::string to_string(const Rational<SafeInt64> &val) {
  int64_t const &val_den = val.get_const_den().get_const_val();
  int64_t const &val_num = val.get_const_num().get_const_val();
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

//
inline void ResInt_Kernel(Rational<SafeInt64> const &a,
                          Rational<SafeInt64> const &b,
                          Rational<SafeInt64> &res) {
  SafeInt64 const &a_den = a.get_const_den();
  SafeInt64 const &b_den = b.get_const_den();
  SafeInt64 const &a_num = a.get_const_num();
  SafeInt64 const &b_num = b.get_const_num();
  SafeInt64 gcd = KernelGcdPair(a_den, b_den);
  SafeInt64 a_mul = a_num * (b_den / gcd);
  SafeInt64 b_mul = b_num * (a_den / gcd);
  SafeInt64 res_i(0);
  ResInt_Kernel(a_mul, b_mul, res_i);
  SafeInt64 prod = a_den * (b_den / gcd);
  res = Rational(res_i, prod);
}

inline SafeInt64 CanonicalizationUnit(SafeInt64 const &eVal) {
  if (eVal < 0)
    return -1;
  return 1;
}

inline Rational<SafeInt64>
CanonicalizationUnit(Rational<SafeInt64> const &eVal) {
  if (eVal < 0)
    return Rational<SafeInt64>(-1);
  return Rational<SafeInt64>(1);
}

inline Rational<SafeInt64> T_NormGen(Rational<SafeInt64> const &x) {
  return T_abs(x);
}

inline SafeInt64 T_NormGen(SafeInt64 const &x) { return T_abs(x); }

inline bool IsInteger(Rational<SafeInt64> const &x) {
  int64_t const &val = x.get_const_den().get_const_val();
  return val == 1;
}

inline Rational<SafeInt64> GetDenominator(Rational<SafeInt64> const &x) {
  SafeInt64 const &eDen = x.get_const_den();
  return Rational(eDen);
}

// We need to have nbRow as input for template reasons. But it is unused in the
// symmetric case. So, pragma statement is needed to avoid a warning being
// thrown.

inline SafeInt64 GetDenominator([[maybe_unused]] SafeInt64 const &x) {
  return 1;
}

inline SafeInt64 GetNumerator_z(Rational<SafeInt64> const &x) {
  return x.get_const_num();
}

inline SafeInt64 GetDenominator_z(Rational<SafeInt64> const &x) {
  return x.get_const_den();
}

inline SafeInt64 GetNumerator_z(SafeInt64 const &x) { return x; }

inline SafeInt64 GetDenominator_z([[maybe_unused]] SafeInt64 const &x) {
  return 1;
}

inline void ScalingInteger_Kernel(stc<Rational<SafeInt64>> const &x,
                                  SafeInt64 &x_ret) {
  x_ret = x.val.get_const_den();
}

inline void ScalingInteger_Kernel([[maybe_unused]] stc<SafeInt64> const &x,
                                  SafeInt64 &x_ret) {
  x_ret = 1;
}

// Rational<SafeInt64> as input

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1,
                            Rational<SafeInt64> &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1, double &a2) {
  int64_t const &val_den = a1.val.get_const_den().get_const_val();
  int64_t const &val_num = a1.val.get_const_num().get_const_val();
  double val_den_d = static_cast<double>(val_den);
  double val_num_d = static_cast<double>(val_num);
  a2 = val_num_d / val_den_d;
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
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  SafeInt64 val2(val1);
  a2 = Rational(val2);
}

inline void TYPE_CONVERSION(stc<uint8_t> const &a1, SafeInt64 &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  a2 = SafeInt64(val1);
}

// uint16_t as input

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, Rational<SafeInt64> &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  SafeInt64 val2(val1);
  a2 = Rational(val2);
}

inline void TYPE_CONVERSION(stc<uint16_t> const &a1, SafeInt64 &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val = static_cast<Tint>(a1.val);
  a2 = SafeInt64(val);
}

// uint32_t as input

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, Rational<SafeInt64> &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  SafeInt64 val2(val1);
  a2 = Rational(val2);
}

inline void TYPE_CONVERSION(stc<uint32_t> const &a1, SafeInt64 &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val = static_cast<Tint>(a1.val);
  a2 = SafeInt64(val);
}

// Conversion to SafeInt64

#ifdef __APPLE__
// long as input (which is not the same as int64_t on APPLE platform)

inline void TYPE_CONVERSION(stc<long> const &a1, Rational<SafeInt64> &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  SafeInt64 val2(val1);
  a2 = Rational(val2);
}

inline void TYPE_CONVERSION(stc<long> const &a1, SafeInt64 &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  a2 = SafeInt64(val1);
}

#endif

// int64_t as input

inline void TYPE_CONVERSION(stc<int64_t> const &a1, Rational<SafeInt64> &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  SafeInt64 val2(val1);
  a2 = Rational(val2);
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, SafeInt64 &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  a2 = SafeInt64(val1);
}

// int32_t as input

inline void TYPE_CONVERSION(stc<int32_t> const &a1, Rational<SafeInt64> &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  SafeInt64 val2(val1);
  a2 = Rational(val2);
}

inline void TYPE_CONVERSION(stc<int32_t> const &a1, SafeInt64 &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  a2 = SafeInt64(val1);
}

// int16_t as input

inline void TYPE_CONVERSION(stc<int16_t> const &a1, Rational<SafeInt64> &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  SafeInt64 val2(val1);
  a2 = Rational(val2);
}

inline void TYPE_CONVERSION(stc<int16_t> const &a1, SafeInt64 &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  a2 = SafeInt64(val1);
}

// int8_t as input

inline void TYPE_CONVERSION(stc<int8_t> const &a1, Rational<SafeInt64> &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  SafeInt64 val2(val1);
  a2 = Rational(val2);
}

inline void TYPE_CONVERSION(stc<int8_t> const &a1, SafeInt64 &a2) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(a1.val);
  a2 = SafeInt64(val1);
}

// SafeInt64 as input

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, SafeInt64 &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, Rational<SafeInt64> &a2) {
  a2 = Rational(a1.val);
}

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, int16_t &a2) {
  int64_t const &eVal = a1.val.get_const_val();
  a2 = static_cast<int16_t>(eVal);
}

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, int32_t &a2) {
  int64_t const &eVal = a1.val.get_const_val();
  a2 = static_cast<int32_t>(eVal);
}

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, int64_t &a2) {
  int64_t const &eVal = a1.val.get_const_val();
  a2 = static_cast<int64_t>(eVal);
}

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, double &a2) {
  a2 = a1.val.get_d();
}

inline void TYPE_CONVERSION(stc<SafeInt64> const &a1, T_uint64_t &a2) {
  int64_t const &eVal = a1.val.get_const_val();
  a2 = static_cast<T_uint64_t>(eVal);
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1, SafeInt64 &a2) {
  T_uint64_t const &val1 = a1.val;
  int64_t val2 = static_cast<int64_t>(val1);
  a2 = SafeInt64(val2);
}

inline void TYPE_CONVERSION(stc<Rational<SafeInt64>> const &a1,
                            T_uint64_t &a2) {
  Termination_rat_safeint_not_integer(a1);
  SafeInt64 a1_z = a1.val.get_const_num();
  a2 = static_cast<T_uint64_t>(a1_z.get_const_val());
}

inline void TYPE_CONVERSION(stc<T_uint64_t> const &a1,
                            Rational<SafeInt64> &a2) {
  T_uint64_t const &val1 = a1.val;
  int64_t val2 = static_cast<int64_t>(val1);
  SafeInt64 val3(val2);
  a2 = Rational(val3);
}

bool universal_square_root(SafeInt64 &ret, SafeInt64 const &val) {
  int64_t const &val_i = val.get_const_val();
  if (val_i < 0) {
    return false;
  }
  double val_d = static_cast<double>(val_i);
  double ret_d = sqrt(val_d);
  int64_t ret_i = static_cast<int64_t>(round(ret_d));
  ret = SafeInt64(ret_i);
  SafeInt64 eProd = ret * ret;
  return eProd == val;
}

bool universal_square_root(Rational<SafeInt64> &ret,
                           Rational<SafeInt64> const &val) {
  SafeInt64 const &val_num = val.get_const_num();
  SafeInt64 const &val_den = val.get_const_den();
  SafeInt64 ret_num, ret_den;
  if (!universal_square_root(ret_num, val_num))
    return false;
  if (!universal_square_root(ret_den, val_den))
    return false;
  ret = Rational(ret_num, ret_den);
  return true;
}

inline void set_to_infinity(Rational<SafeInt64> &x) {
  int64_t val1 = std::numeric_limits<int64_t>::max();
  SafeInt64 val2(val1);
  x = Rational(val2);
}

//
// Nearest integer and similar stuff.
//
inline Rational<SafeInt64> FractionalPart(Rational<SafeInt64> const &x) {
  SafeInt64 const &eNum = x.get_const_num();
  SafeInt64 const &eDen = x.get_const_den();
  SafeInt64 r(0);
  ResInt_Kernel(eNum, eDen, r);
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

inline void FloorInteger(Rational<SafeInt64> const &xI,
                         Rational<SafeInt64> &xO) {
  xO = Floor_safe_rat(xI);
}

inline void FloorInteger(Rational<SafeInt64> const &xI, SafeInt64 &xO) {
  Rational<SafeInt64> xO_q = Floor_safe_rat(xI);
  xO = xO_q.get_const_num();
}

inline void FloorInteger(Rational<SafeInt64> const &xI, int &xO) {
  Rational<SafeInt64> xO_q = Floor_safe_rat(xI);
  xO = static_cast<int>(xO_q.get_const_num().get_const_val());
}

inline void FloorInteger(Rational<SafeInt64> const &xI, long &xO) {
  Rational<SafeInt64> xO_q = Floor_safe_rat(xI);
  xO = static_cast<long>(xO_q.get_const_num().get_const_val());
}

inline void CeilInteger(Rational<SafeInt64> const &xI,
                        Rational<SafeInt64> &xO) {
  xO = Ceil_safe_rat(xI);
}

inline void CeilInteger(Rational<SafeInt64> const &xI, SafeInt64 &xO) {
  Rational<SafeInt64> xO_q = Ceil_safe_rat(xI);
  xO = xO_q.get_const_num();
}

inline void CeilInteger(Rational<SafeInt64> const &xI, int &xO) {
  Rational<SafeInt64> xO_q = Ceil_safe_rat(xI);
  xO = static_cast<int>(xO_q.get_const_num().get_const_val());
}

inline void CeilInteger(Rational<SafeInt64> const &xI, long &xO) {
  Rational<SafeInt64> xO_q = Ceil_safe_rat(xI);
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

inline void NearestInteger(Rational<SafeInt64> const &xI,
                           Rational<SafeInt64> &xO) {
  xO = NearestInteger_rni(xI);
}

inline void NearestInteger(Rational<SafeInt64> const &xI, SafeInt64 &xO) {
  Rational<SafeInt64> xO_q = NearestInteger_rni(xI);
  xO = xO_q.get_const_num();
}

inline void NearestInteger(int const &xI, Rational<SafeInt64> &xO) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(xI);
  SafeInt64 val2(val1);
  xO = Rational(val2);
}

inline void NearestInteger(long const &xI, Rational<SafeInt64> &xO) {
  using Tint = typename SafeInt64::Tint;
  Tint val1 = static_cast<Tint>(xI);
  SafeInt64 val2(val1);
  xO = Rational(val2);
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

// SafeInt64

template <class Archive>
inline void serialize(Archive &ar, SafeInt64 & x,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("the_val", x.get_val());
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYSAFEINT_H_
// clang-format on
