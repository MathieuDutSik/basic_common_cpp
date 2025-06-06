// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_RATIONAL_H_
#define SRC_NUMBER_RATIONAL_H_

// clang-format off
#include "BasicNumberTypes.h"
#include "ExceptionsFunc.h"
#include "ResidueQuotient.h"
#include "TemplateTraits.h"
#include "TypeConversion.h"
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
// clang-format on

/*
  A rational class supposed to build integers.
  Typical use would be for example Rational<long> which should be pretty fast.
 */
template <typename Tint> struct Rational {
private:
  Tint num;
  Tint den;

public:
  // Constructors
  Rational() : num(0), den(1) {}
  Rational(Tint const &x) : num(x), den(1) {}
  Rational(Tint const &num, Tint const &den) : num(num), den(den) {}
  Rational(Rational<Tint> const &x) : num(x.num), den(x.den) {}
  // Assignment operators
  Rational<Tint> operator=(Tint const &u) {
    // assignment operator from int
    num = u;
    den = 1;
    return *this;
  }
  Rational<Tint> operator=(Rational<Tint> const &x) {
    // assignment operator
    num = x.num;
    den = x.den;
    return *this;
  }
  Tint &get_num() { return num; }
  Tint &get_den() { return den; }
  const Tint &get_const_num() const { return num; }
  const Tint &get_const_den() const { return den; }

private:
  // A few internal functions.
  static Tint comp_gcd(Tint const &m, Tint const &n) {
    Tint f = m;
    if (m < 0)
      f = -m;
    Tint g = n;
    if (n < 0)
      g = -n;
    Tint h;
    while (g != 0) {
      h = g;
      g = f % g;
      f = h;
    }
    return f;
  }
  void gcd_reduction() {
    Tint gcd = comp_gcd(num, den);
    num = num / gcd;
    den = den / gcd;
  }
  void check(std::string const &str_test) const {
    if (den < 0) {
      std::cerr << "Error occurred at " << str_test << "\n";
      std::cerr << "num=" << num << " den=" << den << "\n";
      throw TerminalException{1};
    }
  }
  static void check_static(Rational<Tint> const &x,
                           std::string const &str_test) {
    if (x.den < 0) {
      std::cerr << "Error occurred at " << str_test << "\n";
      std::cerr << "x.num=" << x.num << " x.den=" << x.den << "\n";
      throw TerminalException{1};
    }
  }

public:
  void operator/=(Rational<Tint> const &x) {
    if (x.num > 0) {
      num *= x.den;
      den *= x.num;
    } else {
      num *= -x.den;
      den *= -x.num;
    }
    //    check("/=");
    gcd_reduction();
  }
  void operator/=(Tint const &x) {
    if (x > 0) {
      den *= x;
    } else {
      num = -num;
      den *= -x;
    }
    //    check("/= (Tint)");
    gcd_reduction();
  }
  void operator+=(Rational<Tint> const &x) {
    Tint gcd = comp_gcd(den, x.den);
    Tint part_prod = x.den / gcd;
    Tint new_den = den * part_prod;
    num = num * part_prod + x.num * (den / gcd);
    den = new_den;
    // Yes, it is needed: example 1/2 + 1/2
    //    check("+=");
    gcd_reduction();
  }
  void operator-=(Rational<Tint> const &x) {
    Tint gcd = comp_gcd(den, x.den);
    Tint part_prod = x.den / gcd;
    Tint new_den = den * part_prod;
    num = num * part_prod - x.num * (den / gcd);
    den = new_den;
    // Yes, it is needed: example 1/2 - 1/2
    //    check("-=");
    gcd_reduction();
  }
  //
  friend Rational<Tint> operator+(Rational<Tint> const &x,
                                  Rational<Tint> const &y) {
    Rational<Tint> z;
    Tint gcd = Rational<Tint>::comp_gcd(x.den, y.den);
    z.den = y.den * x.den / gcd;
    z.num = y.num * (x.den / gcd) + x.num * (y.den / gcd);
    z.gcd_reduction();
    return z;
  }
  friend Rational<Tint> operator+(int const &x, Rational<Tint> const &y) {
    Rational<Tint> z;
    z.num = x * y.den + y.num;
    z.den = y.den;
    return z;
  }
  friend Rational<Tint> operator+(Rational<Tint> const &x, int const &y) {
    Rational<Tint> z;
    z.num = x.num + y * x.den;
    z.den = x.den;
    return z;
  }
  //
  friend Rational<Tint> operator-(Rational<Tint> const &x,
                                  Rational<Tint> const &y) {
    Rational<Tint> z;
    Tint gcd = Rational<Tint>::comp_gcd(x.den, y.den);
    z.den = x.den * y.den / gcd;
    z.num = x.num * (y.den / gcd) - y.num * (x.den / gcd);
    z.gcd_reduction();
    return z;
  }
  friend Rational<Tint> operator-(int const &x, Rational<Tint> const &y) {
    Rational<Tint> z;
    z.num = x * y.den - y.num;
    z.den = y.den;
    return z;
  }
  friend Rational<Tint> operator-(Rational<Tint> const &x, int const &y) {
    Rational<Tint> z;
    z.num = x.num - y * x.den;
    z.den = x.den;
    return z;
  }
  //
  friend Rational<Tint> operator-(Rational<Tint> const &x) {
    Rational<Tint> z;
    z.num = -x.num;
    z.den = x.den;
    return z;
  }
  friend Rational<Tint> operator/(Tint const &x, Rational<Tint> const &y) {
    Rational<Tint> z;
    if (y.num > 0) {
      z.num = x * y.den;
      z.den = y.num;
    } else {
      z.num = -x * y.den;
      z.den = -y.num;
    }
    z.gcd_reduction();
    return z;
  }
  friend Rational<Tint> operator/(Rational<Tint> const &x,
                                  Rational<Tint> const &y) {
    Rational<Tint> z;
    if (y.num > 0) {
      z.num = x.num * y.den;
      z.den = x.den * y.num;
    } else {
      z.num = -x.num * y.den;
      z.den = -x.den * y.num;
    }
    z.gcd_reduction();
    return z;
  }
  void operator*=(Rational<Tint> const &x) {
    num = num * x.num;
    den = den * x.den;
    //    check("*");
    gcd_reduction();
  }
  //
  friend Rational<Tint> operator*(Rational<Tint> const &x,
                                  Rational<Tint> const &y) {
    Rational<Tint> z;
    z.num = x.num * y.num;
    z.den = x.den * y.den;
    z.gcd_reduction();
    return z;
  }
  friend Rational<Tint> operator*(int const &x, Rational<Tint> const &y) {
    Rational<Tint> z;
    z.num = x * y.num;
    z.den = y.den;
    z.gcd_reduction();
    return z;
  }
  friend Rational<Tint> operator*(Rational<Tint> const &x, int const &y) {
    Rational<Tint> z;
    z.num = x.num * y;
    z.den = x.den;
    z.gcd_reduction();
    return z;
  }
  //
  friend std::ostream &operator<<(std::ostream &os, Rational<Tint> const &v) {
    if (v.den == 1)
      return os << v.num;
    else
      return os << v.num << "/" << v.den;
  }
  friend std::istream &operator>>(std::istream &is, Rational<Tint> &v) {
    char c;
    std::string s;
    size_t miss_val = std::numeric_limits<size_t>::max();
    size_t pos_slash = miss_val;
    size_t pos = 0;
    // First skipping the spaces.
    while (true) {
      // is.get(c) will read characters but is >> c skip the spaces.
      is.get(c);
      if (c != ' ' && c != '\n') {
        s += c;
        // First character cannot be a slash
        pos++;
        break;
      }
    }
    while (true) {
      if (is.eof()) {
        break;
      }
      is.get(c);
      if (c == ' ' || c == '\n') {
        break;
      }
      s += c;
      if (c == '/')
        pos_slash = pos;
      pos++;
    }
    // It seems the istringstream is not working properly in some cases.
    // For example
    // std::istringstream is("3")
    // int test_val
    // is >> test_val
    // which gets us a value of 33 for test_val
    if (pos_slash == miss_val) {
      std::istringstream isN(s);
      isN >> v.num;
      v.den = 1;
    } else {
      std::string sN = s.substr(0, pos_slash);
      std::string sD = s.substr(pos_slash + 1, s.size() - 1 - pos_slash);
      std::istringstream isN(sN);
      isN >> v.num;
      std::istringstream isD(sD);
      isD >> v.den;
    }
    return is;
  }
  //
  friend bool operator==(Rational<Tint> const &x, Rational<Tint> const &y) {
    return (x.num == y.num) && (x.den == y.den);
  }
  friend bool operator==(Rational<Tint> const &x, int const &y) {
    return (x.num == y) && (x.den == 1);
  }
  //
  friend bool operator!=(Rational<Tint> const &x, Rational<Tint> const &y) {
    return (x.num != y.num) || (x.den != y.den);
  }
  friend bool operator!=(Rational<Tint> const &x, int const &y) {
    if (x.den > 1)
      return true;
    return x.num != y;
  }
  //
  friend bool operator>=(Rational<Tint> const &x, Rational<Tint> const &y) {
    // x >= y is equivalent to x_n * y_d >= y_n * x_d
    return x.num * y.den >= y.num * x.den;
  }
  friend bool operator>=(Rational<Tint> const &x, int const &y) {
    return x.num >= y * x.den;
  }
  //
  friend bool operator<=(Rational<Tint> const &x, Rational<Tint> const &y) {
    return x.num * y.den <= y.num * x.den;
  }
  friend bool operator<=(Rational<Tint> const &x, int const &y) {
    return x.num <= y * x.den;
  }
  //
  friend bool operator>(Rational<Tint> const &x, Rational<Tint> const &y) {
    return x.num * y.den > y.num * x.den;
  }
  friend bool operator>(Rational<Tint> const &x, int const &y) {
    return x.num > y * x.den;
  }
  //
  friend bool operator<(Rational<Tint> const &x, Rational<Tint> const &y) {
    return x.num * y.den < y.num * x.den;
  }
  friend bool operator<(Rational<Tint> const &x, Tint const &y) {
    return x.num < y * x.den;
  }
};

namespace std {
template <typename Tint> std::string to_string(const Rational<Tint> &e_val) {
  std::stringstream s;
  s << e_val;
  std::string converted(s.str());
  return converted;
}
// clang-format off
}  // namespace std
// clang-format on

template <typename Tint> struct is_euclidean_domain<Rational<Tint>> {
  static const bool value = true;
};

template <typename Tint> struct is_exact_arithmetic<Rational<Tint>> {
  static const bool value = is_exact_arithmetic<Tint>::value;
};

template <typename Tint> struct is_implementation_of_Z<Rational<Tint>> {
  static const bool value = false;
};

template <typename Tint> struct is_implementation_of_Q<Rational<Tint>> {
  static const bool value = true;
};

template <typename Tint> struct is_ring_field<Rational<Tint>> {
  static const bool value = true;
};

template <typename Tint> struct is_totally_ordered<Rational<Tint>> {
  static const bool value = true;
};

template <typename Tint> struct underlying_ring<Rational<Tint>> {
  typedef Tint ring_type;
};

template <typename Tint>
struct underlying_totally_ordered_ring<Rational<Tint>> {
  typedef Rational<Tint> real_type;
};

template <typename Tint>
inline Rational<Tint> T_NormGen(Rational<Tint> const &x) {
  return T_abs(x);
}

template <typename Tint>
inline std::pair<Rational<Tint>, Rational<Tint>>
ResQuoInt_kernel(Rational<Tint> const &a, Rational<Tint> const &b) {
  // a = a_n / a_d
  // b = b_n / b_d
  // a = res + q * b  with 0 <= res < |b|
  // equivalent to
  // a_n / a_d = res + q * (b_n / b_d)
  // equivalent to
  // a_n * b_d = res * a_d * b_d + (q * a_d) * b_n
  using Tf = Rational<Tint>;
  Tint a_n = a.get_const_num();
  Tint b_n = b.get_const_num();
  Tint a_d = a.get_const_den();
  Tint b_d = b.get_const_den();
  Tint a1 = a_n * b_d;
  Tint b1 = a_d * b_n;
  Tint q = a1 / b1;
  Tf res = a - q * b;
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
          std::cerr << "Error in ResQuoInt_kernel for raional\n";
          std::cerr << "a=" << a << " b=" << b << "\n";
          std::cerr << "res=" << res << " q=" << q << "\n";
          Tf val = res + q * b;
          Tf val1 = q * b;
          Tf val2 = res + val1;
          std::cerr << "v=" << val << " val1=" << val1 << " val2=" << val2
                    << "\n";
          throw TerminalException{1};
        }
        // Conversion of q from Tint to Rational<Tint> occurring here.
        return {res, q};
      }
    }
  }
}

template <typename Tint>
inline void ResInt_Kernel(Rational<Tint> const &a, Rational<Tint> const &b,
                          Rational<Tint> &res) {
  res = ResQuoInt_kernel(a, b).first;
}

template <typename Tint>
inline void QUO_INT(stc<Rational<Tint>> const &a, stc<Rational<Tint>> const &b,
                    Rational<Tint> &q) {
  q = ResQuoInt_kernel(a.val, b.val).second;
}

#include "QuoIntFcts.h"

// Overlying fields, the original motivation

template <> struct overlying_field<int8_t> {
  typedef Rational<int8_t> field_type;
};

template <> struct overlying_field<int16_t> {
  typedef Rational<int16_t> field_type;
};

template <> struct overlying_field<int32_t> {
  typedef Rational<int32_t> field_type;
};

template <> struct overlying_field<int64_t> {
  typedef Rational<int64_t> field_type;
};

template <> struct overlying_field<Rational<int8_t>> {
  typedef Rational<int8_t> field_type;
};

template <> struct overlying_field<Rational<int16_t>> {
  typedef Rational<int16_t> field_type;
};

template <> struct overlying_field<Rational<int32_t>> {
  typedef Rational<int32_t> field_type;
};

template <> struct overlying_field<Rational<int64_t>> {
  typedef Rational<int64_t> field_type;
};

// hashing function

namespace std {
template <typename T> struct hash<Rational<T>> {
  std::size_t operator()(const Rational<T> &x) const {
    auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
      seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };
    size_t seed = std::hash<T>()(x.get_const_num());
    size_t e_hash = std::hash<T>()(x.get_const_den());
    combine_hash(seed, e_hash);
    return seed;
  }
};
// clang-format off
}  // namespace std
// clang-format on

template <typename T>
void TYPE_CONVERSION_Rational_T(stc<Rational<T>> const &a1, T &a2) {
  const T &den = a1.val.get_const_den();
  if (den != 1) {
    std::string str_err =
        "The denominator should be 1. It is den = " + std::to_string(den);
    throw ConversionException{str_err};
  }
  a2 = a1.val.get_const_num();
}

// The conversion tools (int8_t)

inline void TYPE_CONVERSION(stc<Rational<int8_t>> const &a1,
                            Rational<int8_t> &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<Rational<int8_t>> const &a1, int8_t &a2) {
  TYPE_CONVERSION_Rational_T<int8_t>(a1, a2);
}

inline void TYPE_CONVERSION(stc<int8_t> const &a1, Rational<int8_t> &a2) {
  a2 = a1.val;
}

// The conversion tools (int16_t)

inline void TYPE_CONVERSION(stc<Rational<int16_t>> const &a1,
                            Rational<int16_t> &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<Rational<int16_t>> const &a1, int16_t &a2) {
  TYPE_CONVERSION_Rational_T<int16_t>(a1, a2);
}

inline void TYPE_CONVERSION(stc<int16_t> const &a1, Rational<int16_t> &a2) {
  a2 = a1.val;
}

// The conversion tools (int32_t)

inline void TYPE_CONVERSION(stc<Rational<int32_t>> const &a1,
                            Rational<int32_t> &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<Rational<int32_t>> const &a1, int32_t &a2) {
  TYPE_CONVERSION_Rational_T<int32_t>(a1, a2);
}

inline void TYPE_CONVERSION(stc<int32_t> const &a1, Rational<int32_t> &a2) {
  a2 = a1.val;
}

// The conversion tools (int64_t)

inline void TYPE_CONVERSION(stc<Rational<int64_t>> const &a1,
                            Rational<int64_t> &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<Rational<int64_t>> const &a1, int64_t &a2) {
  TYPE_CONVERSION_Rational_T<int64_t>(a1, a2);
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, Rational<int64_t> &a2) {
  a2 = a1.val;
}

// Obtention of denominators / numerators

template <typename Tint>
inline Rational<Tint> GetNumerator(Rational<Tint> const &x) {
  return x.get_const_num();
}

template <typename Tint>
inline Rational<Tint> GetDenominator(Rational<Tint> const &x) {
  return x.get_const_den();
}

template <typename Tint> inline Tint GetNumerator_z(Rational<Tint> const &x) {
  return x.get_const_num();
}

template <typename Tint> inline Tint GetDenominator_z(Rational<Tint> const &x) {
  return x.get_const_den();
}

template <typename Tint>
void ScalingInteger_Kernel(stc<Rational<Tint>> const &x, Tint &x_ret) {
  x_ret = x.val.get_const_den();
}

// Floor / Ceil / Nearest operations

template <typename Tint>
Rational<Tint> FractionalPart(Rational<Tint> const &x) {
  Tint res = ResInt(x.get_const_num(), x.get_const_den());
  Rational<Tint> fr(res, x.get_const_den());
  return fr;
}

template <typename Tint>
inline void FloorInteger(Rational<Tint> const &xI, Rational<Tint> &xO) {
  Rational<Tint> fr = FractionalPart(xI);
  xO = xI - fr;
}

template <typename Tint>
inline void CeilInteger(Rational<Tint> const &xI, Rational<Tint> &xO) {
  Rational<Tint> fr = FractionalPart(xI);
  if (fr == 0) {
    xO = xI;
  } else {
    xO = 1 + xI - fr;
  }
}

template <typename Tint>
Rational<Tint> NearestInteger_rni(Rational<Tint> const &a) {
  Rational<Tint> fr = FractionalPart(a);
  Rational<Tint> eDiff1 = fr;
  Rational<Tint> eDiff2 = 1 - fr;
  Rational<Tint> RetVal = a - fr;
  if (eDiff1 <= eDiff2) {
    return RetVal;
  } else {
    return 1 + RetVal;
  }
}

template <typename Tint>
inline void NearestInteger(Rational<Tint> const &xI, Rational<Tint> &xO) {
  xO = NearestInteger_rni(xI);
}

template <typename Tint>
inline void NearestInteger(Rational<Tint> const &xI, Tint &xO) {
  Rational<Tint> xO_q = NearestInteger_rni(xI);
  xO = xO_q.get_num();
}

// Serialization stuff

namespace boost::serialization {

template <class Archive, typename T>
inline void serialize(Archive &ar, Rational<T> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("rational_num", val.get_num());
  ar &make_nvp("rational_den", val.get_den());
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

// clang-format off
#endif  // SRC_NUMBER_RATIONAL_H_
// clang-format on
