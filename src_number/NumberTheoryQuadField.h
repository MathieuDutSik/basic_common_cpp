// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYQUADFIELD_H_
#define SRC_NUMBER_NUMBERTHEORYQUADFIELD_H_
// clang-format off
#include "Temp_common.h"
#include "InputOutput.h"
#include <limits>
#include <string>
// clang-format on

template <typename Tinp, int d> class QuadField {
public:
  using Tresidual = Tinp;

private:
  using T = Tinp;
  T a;
  T b;

public:
  T &get_a() { return a; }
  T &get_b() { return b; }
  const T &get_const_a() const { return a; }
  const T &get_const_b() const { return b; }

  // Note: We are putting "int" as argument here because we want to do the
  // comparison with the stuff like x > 0 or x = 1. For the type "rational<T>"
  // we had to forbid that because this lead to erroneous conversion of say
  // int64_t to int with catastrophic loss of precision. But for the
  // QuadField<T> the loss of precision does not occur because T is typically
  // mpq_class. or some other type that does not convert to integers easily. And
  // at the same time the natural conversion of int to int64_t allows the
  // comparison x > 0 and equality set x = 1 to work despite the lack of a
  // operator=(int const& u)

  // Constructor
  QuadField() : a(0), b(0) {}
  QuadField(int const &u) : a(u), b(0) {}
  QuadField(T const &u) : a(u), b(0) {}
  QuadField(T const &_a, T const &_b) : a(_a), b(_b) {}
  QuadField(QuadField<T, d> const &x) : a(x.a), b(x.b) {}
  //  QuadField<T,d>& operator=(QuadField<T,d> const&); // assignment operator
  //  QuadField<T,d>& operator=(T const&); // assignment operator from T
  //  QuadField<T,d>& operator=(int const&); // assignment operator from T
  // assignment operator from int
  QuadField<T, d> operator=(int const &u) {
    a = u;
    b = 0;
    return *this;
  }
  // assignment operator
  QuadField<T, d> operator=(QuadField<T, d> const &x) {
    a = x.a;
    b = x.b;
    return *this;
  }
  //
  // Arithmetic operators below:
  void operator+=(QuadField<T, d> const &x) {
    a += x.a;
    b += x.b;
  }
  void operator-=(QuadField<T, d> const &x) {
    a -= x.a;
    b -= x.b;
  }
  void operator/=(QuadField<T, d> const &x) {
    T disc = x.a * x.a - d * x.b * x.b;
    T a_new = (a * x.a - d * b * x.b) / disc;
    b = (b * x.a - a * x.b) / disc;
    a = a_new;
  }
  friend QuadField<T, d> operator+(QuadField<T, d> const &x,
                                   QuadField<T, d> const &y) {
    return QuadField<T, d>(x.a + y.a, x.b + y.b);
  }
  friend QuadField<T, d> operator-(QuadField<T, d> const &x,
                                   QuadField<T, d> const &y) {
    return QuadField<T, d>(x.a - y.a, x.b - y.b);
  }
  friend QuadField<T, d> operator-(QuadField<T, d> const &x, int const &y) {
    return QuadField<T, d>(x.a - y, x.b);
  }
  friend QuadField<T, d> operator-(QuadField<T, d> const &x) {
    return QuadField<T, d>(-x.a, -x.b);
  }
  friend QuadField<T, d> operator/(int const &x, QuadField<T, d> const &y) {
    QuadField<T, d> z;
    T disc = y.a * y.a - d * y.b * y.b;
    z.a = x * y.a / disc;
    z.b = -x * y.b / disc;
    return z;
  }
  friend QuadField<T, d> operator/(QuadField<T, d> const &x,
                                   QuadField<T, d> const &y) {
    QuadField<T, d> z;
    T disc = y.a * y.a - d * y.b * y.b;
    z.a = (x.a * y.a - d * x.b * y.b) / disc;
    z.b = (x.b * y.a - x.a * y.b) / disc;
    return z;
  }
  void operator*=(QuadField<T, d> const &x) {
    T hA = a * x.a + d * b * x.b;
    b = a * x.b + b * x.a;
    a = hA;
  }
  friend QuadField<T, d> operator*(QuadField<T, d> const &x,
                                   QuadField<T, d> const &y) {
    return QuadField<T, d>(x.a * y.a + d * x.b * y.b, x.a * y.b + x.b * y.a);
  }
  friend QuadField<T, d> operator*(int const &x, QuadField<T, d> const &y) {
    return QuadField<T, d>(x * y.a, x * y.b);
  }
  friend std::ostream &operator<<(std::ostream &os, QuadField<T, d> const &v) {
    std::vector<T> V{v.a, v.b};
    WriteVectorFromRealAlgebraicString(os, V);
    return os;
  }
  friend std::istream &operator>>(std::istream &is, QuadField<T, d> &v) {
    std::vector<T> V = ReadVectorFromRealAlgebraicString<T>(is, 2);
    v.a = V[0];
    v.b = V[1];
    return is;
  }
  friend bool operator==(QuadField<T, d> const &x, QuadField<T, d> const &y) {
    if (x.a != y.a)
      return false;
    if (x.b != y.b)
      return false;
    return true;
  }
  friend bool operator!=(QuadField<T, d> const &x, QuadField<T, d> const &y) {
    if (x.a != y.a)
      return true;
    if (x.b != y.b)
      return true;
    return false;
  }
  friend bool operator!=(QuadField<T, d> const &x, int const &y) {
    if (x.a != y)
      return true;
    if (x.b != 0)
      return true;
    return false;
  }
  friend bool IsNonNegative(QuadField<T, d> const &x) {
    if (x.a == 0 && x.b == 0)
      return true;
    if (x.a >= 0 && x.b >= 0)
      return true;
    if (x.a <= 0 && x.b <= 0)
      return false;
    T disc = x.a * x.a - d * x.b * x.b;
    if (disc > 0) {
      if (x.a >= 0 && x.b <= 0)
        return true;
      if (x.a <= 0 && x.b >= 0)
        return false;
    } else {
      if (x.a >= 0 && x.b <= 0)
        return false;
      if (x.a <= 0 && x.b >= 0)
        return true;
    }
    std::cerr << "Major errors in the code\n";
    return false;
  }
  friend bool operator>=(QuadField<T, d> const &x, QuadField<T, d> const &y) {
    QuadField<T, d> z;
    z = x - y;
    return IsNonNegative(z);
  }
  friend bool operator>=(QuadField<T, d> const &x, int const &y) {
    QuadField<T, d> z;
    z = x - y;
    return IsNonNegative(z);
  }
  friend bool operator<=(QuadField<T, d> const &x, QuadField<T, d> const &y) {
    QuadField<T, d> z;
    z = y - x;
    return IsNonNegative(z);
  }
  friend bool operator<=(QuadField<T, d> const &x, int const &y) {
    QuadField<T, d> z;
    z = y - x;
    return IsNonNegative(z);
  }
  friend bool operator>(QuadField<T, d> const &x, QuadField<T, d> const &y) {
    QuadField<T, d> z;
    z = x - y;
    if (z.a == 0 && z.b == 0)
      return false;
    return IsNonNegative(z);
  }
  friend bool operator>(QuadField<T, d> const &x, int const &y) {
    QuadField<T, d> z;
    z = x - y;
    if (z.a == 0 && z.b == 0)
      return false;
    return IsNonNegative(z);
  }
  friend bool operator<(QuadField<T, d> const &x, QuadField<T, d> const &y) {
    QuadField<T, d> z;
    z = y - x;
    if (z.a == 0 && z.b == 0)
      return false;
    return IsNonNegative(z);
  }
  friend bool operator<(QuadField<T, d> const &x, int const &y) {
    QuadField<T, d> z;
    z = y - x;
    if (z.a == 0 && z.b == 0)
      return false;
    return IsNonNegative(z);
  }
};

template <typename T, int d> struct overlying_field<QuadField<T, d>> {
  typedef QuadField<typename QuadField<T, d>::Tresidual, d> field_type;
};

// Note that the underlying ring is not unique, there are many possibiliies
// actually but we can represent only one in our scheme.
template <typename T, int d> struct underlying_ring<QuadField<T, d>> {
  typedef QuadField<typename QuadField<T, d>::Tresidual, d> ring_type;
};

template <typename T, int d>
inline void TYPE_CONVERSION(stc<QuadField<T, d>> const &x1, double &x2) {
  stc<T> a1{x1.val.get_const_a()};
  stc<T> b1{x1.val.get_const_b()};
  double a2, b2;
  TYPE_CONVERSION(a1, a2);
  TYPE_CONVERSION(b1, b2);
  x2 = a2 + sqrt(d) * b2;
}

template <typename T, int d> struct is_totally_ordered<QuadField<T, d>> {
  static const bool value = true;
};

template <typename T, int d> struct is_ring_field<QuadField<T, d>> {
  static const bool value = is_ring_field<T>::value;
};

template <typename T, int d> struct is_exact_arithmetic<QuadField<T, d>> {
  static const bool value = true;
};

// Hashing function

template <typename T, int d> struct is_implementation_of_Z<QuadField<T, d>> {
  static const bool value = false;
};

template <typename T, int d> struct is_implementation_of_Q<QuadField<T, d>> {
  static const bool value = false;
};

// Hashing function

namespace std {
template <typename T, int d> struct hash<QuadField<T, d>> {
  std::size_t operator()(const QuadField<T, d> &x) const {
    auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
      seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };
    size_t seed = std::hash<int>()(d);
    size_t e_hash1 = std::hash<T>()(x.get_const_a());
    size_t e_hash2 = std::hash<T>()(x.get_const_b());
    combine_hash(seed, e_hash1);
    combine_hash(seed, e_hash2);
    return seed;
  }
};
// clang-format off
}  // namespace std
// clang-format on

// Local typing info

template <typename T> struct is_quad_field {};

template <typename T, int d> struct is_quad_field<QuadField<T, d>> {
  static const bool value = true;
};

// Some functionality

template <typename T, int d> bool IsInteger(QuadField<T, d> const &x) {
  if (x.get_const_b() != 0)
    return false;
  return IsInteger(x.get_const_a());
}

// The conversion tools (int)

template <typename T1, typename T2, int d>
inline void TYPE_CONVERSION(stc<QuadField<T1, d>> const &x1,
                            QuadField<T2, d> &x2) {
  stc<T1> a1{x1.val.get_const_a()};
  stc<T1> b1{x1.val.get_const_b()};
  TYPE_CONVERSION(a1, x2.get_a());
  TYPE_CONVERSION(b1, x2.get_b());
}

template <typename T1, typename T2, int d>
inline typename std::enable_if<not is_quad_field<T2>::value, void>::type
TYPE_CONVERSION(stc<QuadField<T1, d>> const &x1, T2 &x2) {
  if (x1.val.b != 0) {
    std::string str = "Conversion error for quadratic field";
    throw ConversionException{str};
  }
  stc<T1> a1{x1.val.get_const_a()};
  TYPE_CONVERSION(a1, x2);
}

// Serialization stuff

namespace boost::serialization {

template <class Archive, typename T, int d>
inline void serialize(Archive &ar, QuadField<T, d> &val,
                      [[maybe_unused]] const unsigned int version) {
  ar &make_nvp("quadfield_a", val.get_a());
  ar &make_nvp("quadfield_b", val.get_b());
}

// clang-format off
}  // namespace boost::serialization
// clang-format on

// Turning into something rational

template <typename Tring, typename T, int d>
void ScalingInteger_Kernel(stc<QuadField<T, d>> const &x, Tring &x_res) {
  using Tfield = T;
  Tfield const &a = x.val.get_const_a();
  Tfield const &b = x.val.get_const_b();
  x_res = LCMpair(GetDenominator_z(a), GetDenominator_z(b));
}

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYQUADFIELD_H_
// clang-format on
