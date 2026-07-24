// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYQUADFIELD_H_
#define SRC_NUMBER_NUMBERTHEORYQUADFIELD_H_
// clang-format off
#include "Temp_common.h"
#include "InputOutput.h"
#include <limits>
#include <string>
// clang-format on

template <typename Tinp, int d> class QuadField;

// Lazy product of two QuadField elements -- a minimal expression template (see
// the analogous RatProd in rational.h). `a * b` returns this proxy; the fast
// sinks evaluate it directly into their own buffers,
//   prod  = a * b;   -> QuadField::operator=(QuadProd)    (in place)
//   acc  += a * b;   -> QuadField::operator+=(QuadProd)   (fused, no wrapper)
//   acc  -= a * b;   -> QuadField::operator-=(QuadProd)   (fused, no wrapper)
//   QuadField r=a*b; -> QuadField(QuadProd)               (fresh, as before)
// and every other use materializes it into a QuadField through the operators
// defined after the class, so results are identical to the eager version. Like
// gmpxx expression templates it holds references: consume within the same
// full-expression, do not bind with `auto` and reuse later.
template <typename Tinp, int d> struct QuadProd {
  QuadField<Tinp, d> const &x;
  QuadField<Tinp, d> const &y;
};

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
  // Construct from a lazy product a*b: a fresh object, as the eager operator*
  // did. Also the implicit QuadProd -> QuadField conversion for every non-sink
  // use.
  QuadField(QuadProd<Tinp, d> const &e)
      : a(e.x.a * e.y.a + d * e.x.b * e.y.b),
        b(e.x.a * e.y.b + e.x.b * e.y.a) {}
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
  // Assign from a lazy product a*b: multiply in place, reusing this->a / this->b
  // (only one temporary, as in operator*=). Aliasing-safe when this == x or y:
  // the new a is computed into a temporary and stored last, and the new b reads
  // this->a before it is overwritten.
  QuadField<T, d> &operator=(QuadProd<Tinp, d> const &e) {
    T na = e.x.a * e.y.a + d * e.x.b * e.y.b;
    b = e.x.a * e.y.b + e.x.b * e.y.a;
    a = na;
    return *this;
  }
  //
  // Arithmetic operators below:
  void operator+=(QuadField<T, d> const &x) {
    a += x.a;
    b += x.b;
  }
  // Fused accumulate of a lazy product: this += a*b, without the wrapper
  // QuadField the eager operator* would build. Product components are formed
  // first (aliasing-safe).
  void operator+=(QuadProd<Tinp, d> const &e) {
    T pa = e.x.a * e.y.a + d * e.x.b * e.y.b;
    T pb = e.x.a * e.y.b + e.x.b * e.y.a;
    a += pa;
    b += pb;
  }
  void operator-=(QuadField<T, d> const &x) {
    a -= x.a;
    b -= x.b;
  }
  // Fused subtract of a lazy product: this -= a*b.
  void operator-=(QuadProd<Tinp, d> const &e) {
    T pa = e.x.a * e.y.a + d * e.x.b * e.y.b;
    T pb = e.x.a * e.y.b + e.x.b * e.y.a;
    a -= pa;
    b -= pb;
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
  // Lazy: returns a QuadProd proxy (see above), evaluated in place by the
  // consumer. Mixed int*QuadField stays eager below.
  friend QuadProd<Tinp, d> operator*(QuadField<T, d> const &x,
                                     QuadField<T, d> const &y) {
    return QuadProd<Tinp, d>{x, y};
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

// ---------------------------------------------------------------------------
// QuadProd (the lazy a*b proxy) as a first-class value. Every use other than the
// in-place sinks above materializes the proxy into a QuadField and delegates to
// the ordinary QuadField operators, so results are identical to the eager
// implementation. Arithmetic operators return QuadField explicitly so that a
// QuadProd produced on the right-hand side is materialized before the operand
// temporaries die.
// ---------------------------------------------------------------------------
template <typename Tinp, int d>
inline QuadField<Tinp, d> const &quad_eval(QuadField<Tinp, d> const &x) {
  return x;
}
template <typename Tinp, int d>
inline QuadField<Tinp, d> quad_eval(QuadProd<Tinp, d> const &e) {
  return QuadField<Tinp, d>(e);
}

#define QUADFIELD_QUADPROD_ARITH(OP)                                           \
  template <typename Tinp, int d>                                              \
  inline QuadField<Tinp, d> operator OP(QuadProd<Tinp, d> const &a,            \
                                        QuadProd<Tinp, d> const &b) {          \
    return quad_eval(a) OP quad_eval(b);                                       \
  }                                                                            \
  template <typename Tinp, int d>                                              \
  inline QuadField<Tinp, d> operator OP(QuadProd<Tinp, d> const &a,            \
                                        QuadField<Tinp, d> const &b) {         \
    return quad_eval(a) OP b;                                                  \
  }                                                                            \
  template <typename Tinp, int d>                                              \
  inline QuadField<Tinp, d> operator OP(QuadField<Tinp, d> const &a,           \
                                        QuadProd<Tinp, d> const &b) {          \
    return a OP quad_eval(b);                                                  \
  }
QUADFIELD_QUADPROD_ARITH(+)
QUADFIELD_QUADPROD_ARITH(-)
QUADFIELD_QUADPROD_ARITH(*)
QUADFIELD_QUADPROD_ARITH(/)
#undef QUADFIELD_QUADPROD_ARITH

#define QUADFIELD_QUADPROD_CMP(OP)                                             \
  template <typename Tinp, int d>                                              \
  inline bool operator OP(QuadProd<Tinp, d> const &a,                          \
                          QuadProd<Tinp, d> const &b) {                        \
    return quad_eval(a) OP quad_eval(b);                                       \
  }                                                                            \
  template <typename Tinp, int d>                                              \
  inline bool operator OP(QuadProd<Tinp, d> const &a,                          \
                          QuadField<Tinp, d> const &b) {                       \
    return quad_eval(a) OP b;                                                  \
  }                                                                            \
  template <typename Tinp, int d>                                              \
  inline bool operator OP(QuadField<Tinp, d> const &a,                         \
                          QuadProd<Tinp, d> const &b) {                        \
    return a OP quad_eval(b);                                                  \
  }                                                                            \
  template <typename Tinp, int d>                                              \
  inline bool operator OP(QuadProd<Tinp, d> const &a, int const &b) {          \
    return quad_eval(a) OP b;                                                  \
  }
QUADFIELD_QUADPROD_CMP(==)
QUADFIELD_QUADPROD_CMP(!=)
QUADFIELD_QUADPROD_CMP(<)
QUADFIELD_QUADPROD_CMP(>)
QUADFIELD_QUADPROD_CMP(<=)
QUADFIELD_QUADPROD_CMP(>=)
#undef QUADFIELD_QUADPROD_CMP

template <typename Tinp, int d>
inline QuadField<Tinp, d> operator-(QuadProd<Tinp, d> const &e) {
  return -quad_eval(e);
}
template <typename Tinp, int d>
inline bool IsNonNegative(QuadProd<Tinp, d> const &e) {
  return IsNonNegative(QuadField<Tinp, d>(e));
}
template <typename Tinp, int d>
inline std::ostream &operator<<(std::ostream &os, QuadProd<Tinp, d> const &e) {
  return os << QuadField<Tinp, d>(e);
}

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

template <typename T, int d>
inline void NearestInteger(QuadField<T, d> const &xI, QuadField<T, d> &xO) {
  using Teff = QuadField<T, d>;
  xO = 0;
  while (true) {
    Teff err = T_abs(xI - xO);
    auto get_move = [&]() -> Teff {
      if (xI > xO) {
        return Teff(1);
      } else {
        return Teff(-1);
      }
    };
    Teff delta = get_move();
    Teff xB = xO + delta;
    Teff err_B = T_abs(xB - xI);
    if (err_B >= err) {
      return;
    }
    xO = xB;
  }
}

template <typename T, int d> struct is_totally_ordered<QuadField<T, d>> {
  static const bool value = true;
};

template <typename T, int d> struct is_ring_field<QuadField<T, d>> {
  static const bool value = is_ring_field<T>::value;
};

// FMA form (see is_fma_prefered). The reused-scratch form is fastest for
// QuadField (measured): operator=(QuadProd) uses one temporary vs
// operator+=(QuadProd)'s two.
template <typename T, int d> struct is_fma_prefered<QuadField<T, d>> {
  static const bool value = false;
};

template <typename T, int d> struct is_exact_arithmetic<QuadField<T, d>> {
  static const bool value = true;
};

// Hashing function

template <typename T, int d> struct is_implementation_of_Z<QuadField<T, d>> {
  static const bool value = false;
};

// A quadratic field inherits Bareiss-eligibility from its base: exact over an
// exact base (e.g. QuadField<mpq_class,d>, where Bareiss wins ~2-4x), and off
// over a floating-point base where numerical pivoting must be preserved.
template <typename T, int d>
struct use_bareiss_for_determinants<QuadField<T, d>> {
  static const bool value = use_bareiss_for_determinants<T>::value;
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
requires (!is_quad_field<T2>::value)
inline void TYPE_CONVERSION(stc<QuadField<T1, d>> const &x1, T2 &x2) {
  if (x1.val.get_const_b() != 0) {
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
