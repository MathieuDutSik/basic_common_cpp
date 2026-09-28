// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYQUADFIELD_H_
#define SRC_NUMBER_NUMBERTHEORYQUADFIELD_H_
// clang-format off
#include "Temp_common.h"
#include "InputOutput.h"
#include "IntegerSearch.h"
#include "MatrixTypes.h"
#include <boost/serialization/nvp.hpp>
#include <cmath>
#include <cstdint>
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
    static thread_local T na, ph;
    na = e.x.a * e.y.a;
    ph = e.x.b * e.y.b;
    ph *= d;
    na += ph;
    ph = e.x.a * e.y.b;
    b = ph;
    ph = e.x.b * e.y.a;
    b += ph;
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
    static thread_local T pa, pb, ph;
    pa = e.x.a * e.y.a;
    ph = e.x.b * e.y.b;
    ph *= d;
    pa += ph;
    pb = e.x.a * e.y.b;
    ph = e.x.b * e.y.a;
    pb += ph;
    a += pa;
    b += pb;
  }
  void operator-=(QuadField<T, d> const &x) {
    a -= x.a;
    b -= x.b;
  }
  // Fused subtract of a lazy product: this -= a*b.
  void operator-=(QuadProd<Tinp, d> const &e) {
    static thread_local T pa, pb, ph;
    pa = e.x.a * e.y.a;
    ph = e.x.b * e.y.b;
    ph *= d;
    pa += ph;
    pb = e.x.a * e.y.b;
    ph = e.x.b * e.y.a;
    pb += ph;
    a -= pa;
    b -= pb;
  }
  // Every division here is x / y = x * conj(y) / N(y), with the norm
  // N(y) = c^2 - d e^2 a scalar of the base type: the coordinates of
  // x * conj(y) are divided by it. Over a field that is exact by
  // construction. Over a ring -- the case of Z[sqrt(d)], which is what
  // underlying_ring returns for a quadratic field -- it need not be: the
  // quotient lies in the ring exactly when the norm divides both
  // coordinates, and there is nothing to fall back on when it does not, so
  // the division reports and throws rather than truncate silently.
  static T DivideByNorm(T const &num, T const &norm) {
    if constexpr (is_ring_field<T>::value) {
      return num / norm;
    } else {
      if (norm == 0) {
        std::cerr << "QUADFIELD: division by zero in Z[sqrt(" << d << ")]\n";
        throw TerminalException{1};
      }
      T quot = num / norm;
      if (quot * norm != num) {
        std::cerr << "QUADFIELD: the quotient is not an element of the ring "
                     "Z[sqrt(" << d << ")]\n";
        std::cerr << "QUADFIELD: the norm " << norm
                  << " does not divide the coordinate " << num << "\n";
        std::cerr << "QUADFIELD: use the overlying field for this quotient\n";
        throw TerminalException{1};
      }
      return quot;
    }
  }
  void operator/=(QuadField<T, d> const &x) {
    T disc = x.a * x.a - d * x.b * x.b;
    T a_new = DivideByNorm(a * x.a - d * b * x.b, disc);
    b = DivideByNorm(b * x.a - a * x.b, disc);
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
    z.a = DivideByNorm(x * y.a, disc);
    z.b = DivideByNorm(-x * y.b, disc);
    return z;
  }
  friend QuadField<T, d> operator/(QuadField<T, d> const &x,
                                   QuadField<T, d> const &y) {
    QuadField<T, d> z;
    T disc = y.a * y.a - d * y.b * y.b;
    z.a = DivideByNorm(x.a * y.a - d * x.b * y.b, disc);
    z.b = DivideByNorm(x.b * y.a - x.a * y.b, disc);
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

// The field of fractions and the underlying ring follow the base type: over
// mpz_class the quadratic field is the order Z[sqrt(d)] and its field of
// fractions is the quadratic field over mpq_class. Both are fixpoints on the
// side they already are, so overlying_field of a field and underlying_ring of
// a ring return the type itself.
template <typename T, int d> struct overlying_field<QuadField<T, d>> {
  typedef QuadField<typename overlying_field<T>::field_type, d> field_type;
};

// The underlying ring is Z[sqrt(d)] = { a + b sqrt(d) : a, b in Z }, the free
// Z-module on 1 and sqrt(d). It is not the ring of integers of the field: for
// d = 1 mod 4 that is the strictly larger Z[(1+sqrt(d))/2], which the (a, b)
// layout over (1, sqrt(d)) cannot even represent. It is not canonical either,
// but it is a ring, sqrt(d) being a root of the monic X^2 - d, and running
// over it avoids the denominators of the field.
template <typename T, int d> struct underlying_ring<QuadField<T, d>> {
  typedef QuadField<typename underlying_ring<T>::ring_type, d> ring_type;
};

// The rational scalars, where underlying_ring above and underlying_z_ring
// part ways: the ring for the fraction-free paths keeps sqrt(d) and is
// Z[sqrt(d)], while the rational integers of Q(sqrt(d)) are plain Z. Code
// working on a lattice over Z -- the basis transformation of an LLL
// reduction, the index of a sublattice being factored -- wants this one, and
// running it over Z[sqrt(d)] would ask a euclidean division of a ring that
// has none for most d.
template <typename T, int d> struct underlying_z_ring<QuadField<T, d>> {
  typedef typename underlying_z_ring<T>::ring_type ring_type;
};

// Only over a base that is a field, Z[sqrt(d)] containing no Q. The requires
// clause is what makes that a plain absence -- the specialization does not
// apply and the empty primary is chosen -- rather than a member typedef that
// is an error to instantiate, which no detection idiom could see.
template <typename T, int d>
requires requires { typename underlying_q_field<T>::field_type; }
struct underlying_q_field<QuadField<T, d>> {
  typedef typename underlying_q_field<T>::field_type field_type;
};

// Answers "is this type a QuadField". The primary has to carry a false value
// rather than be left empty: the trait is read from the requires clause of the
// conversion out of QuadField below, and on a type that never specializes it a
// missing value member is a substitution failure, which drops that conversion
// from the overload set instead of selecting it.
template <typename T> struct is_quad_field {
  static const bool value = false;
};

template <typename T, int d> struct is_quad_field<QuadField<T, d>> {
  static const bool value = true;
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

// See T_frexp in TypeConversion.h. For x = a + b sqrt(d) with a and b of the
// same sign the two terms are added. With opposite signs they would cancel,
// so the exact identity x = (a^2 - d b^2) / (a - b sqrt(d)) is used instead:
// the norm a^2 - d b^2 is computed exactly in T, and the denominator adds two
// terms of the same sign. Either way no digit is lost to cancellation and the
// result has the accuracy of a double.
template <typename T, int d>
inline double T_frexp(QuadField<T, d> const &x, long &e) {
  T const &a = x.get_const_a();
  T const &b = x.get_const_b();
  long e_a, e_b;
  double m_a = T_frexp(a, e_a);
  double m_b = T_frexp(b, e_b) * std::sqrt(static_cast<double>(d));
  bool same_sign = (a >= 0 && b >= 0) || (a <= 0 && b <= 0);
  if (same_sign) {
    return T_frexp_sum(m_a, e_a, m_b, e_b, e);
  }
  T norm = a * a - T(d) * b * b;
  long e_norm, e_conj;
  double m_norm = T_frexp(norm, e_norm);
  double m_conj = T_frexp_sum(m_a, e_a, -m_b, e_b, e_conj);
  return T_frexp_quotient(m_norm, e_norm, m_conj, e_conj, e);
}

// The rounding of an element of Q(sqrt(d)) to an element of Z, in the three
// output types the callers ask for: the field itself, the ring Z[sqrt(d)] that
// underlying_ring makes its ring, and a plain integer type. The lattice an LLL
// reduction works on is Z^n whatever field the form takes its values in, so
// the basis transformation stays in Z and the reduction asks for the last of
// the three.
//
// Each rounding is one exact search of helper_largest_integer_satisfying,
// comparing xI to integers, with no approximation of sqrt(d):
//   --- the floor, the largest n in Z with n <= xI,
//   --- the ceiling, the smallest n in Z with xI <= n,
//   --- the nearest integer, a tie y + 1/2 going to y, as for mpq_class.
// All three commute with the integer translations, Floor(x + n) = Floor(x) +
// n and so on. The helpers return the integer; the overloads below put it in
// the output type.
template <typename T, int d>
inline double helper_quad_field_double(QuadField<T, d> const &xI) {
  double x_d;
  TYPE_CONVERSION(stc<QuadField<T, d>>{xI}, x_d);
  return x_d;
}

template <typename T, int d>
inline typename underlying_z_ring<T>::ring_type
helper_quad_field_floor(QuadField<T, d> const &xI) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  auto pred = [&](Tint const &n) -> bool {
    return QuadField<T, d>(UniversalScalarConversion<T, Tint>(n)) <= xI;
  };
  return helper_largest_integer_satisfying<Tint>(
      std::floor(helper_quad_field_double(xI)), pred);
}

template <typename T, int d>
inline typename underlying_z_ring<T>::ring_type
helper_quad_field_ceil(QuadField<T, d> const &xI) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  // The largest n with n < xI, one below the ceiling.
  auto pred = [&](Tint const &n) -> bool {
    return QuadField<T, d>(UniversalScalarConversion<T, Tint>(n)) < xI;
  };
  Tint n = helper_largest_integer_satisfying<Tint>(
      std::ceil(helper_quad_field_double(xI)) - 1, pred);
  return n + 1;
}

template <typename T, int d>
inline typename underlying_z_ring<T>::ring_type
helper_quad_field_nearest(QuadField<T, d> const &xI) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  // The largest n with n - 1/2 < xI, that is 2 n - 1 < 2 xI: T may be a ring
  // with no 1/2 in it.
  QuadField<T, d> x2 = xI + xI;
  auto pred = [&](Tint const &n) -> bool {
    Tint num = 2 * n - 1;
    return QuadField<T, d>(UniversalScalarConversion<T, Tint>(num)) < x2;
  };
  return helper_largest_integer_satisfying<Tint>(
      std::floor(helper_quad_field_double(xI) + 0.5), pred);
}

template <typename T, int d>
inline void FloorInteger(QuadField<T, d> const &xI, QuadField<T, d> &xO) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  xO = QuadField<T, d>(
      UniversalScalarConversion<T, Tint>(helper_quad_field_floor(xI)));
}

template <typename T, typename Tring, int d>
inline void FloorInteger(QuadField<T, d> const &xI, QuadField<Tring, d> &xO) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  xO = QuadField<Tring, d>(
      UniversalScalarConversion<Tring, Tint>(helper_quad_field_floor(xI)));
}

template <typename T, typename Tout, int d>
requires (!is_quad_field<Tout>::value)
inline void FloorInteger(QuadField<T, d> const &xI, Tout &xO) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  xO = UniversalScalarConversion<Tout, Tint>(helper_quad_field_floor(xI));
}

template <typename T, int d>
inline void CeilInteger(QuadField<T, d> const &xI, QuadField<T, d> &xO) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  xO = QuadField<T, d>(
      UniversalScalarConversion<T, Tint>(helper_quad_field_ceil(xI)));
}

template <typename T, typename Tring, int d>
inline void CeilInteger(QuadField<T, d> const &xI, QuadField<Tring, d> &xO) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  xO = QuadField<Tring, d>(
      UniversalScalarConversion<Tring, Tint>(helper_quad_field_ceil(xI)));
}

template <typename T, typename Tout, int d>
requires (!is_quad_field<Tout>::value)
inline void CeilInteger(QuadField<T, d> const &xI, Tout &xO) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  xO = UniversalScalarConversion<Tout, Tint>(helper_quad_field_ceil(xI));
}

template <typename T, int d>
inline void NearestInteger(QuadField<T, d> const &xI, QuadField<T, d> &xO) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  xO = QuadField<T, d>(
      UniversalScalarConversion<T, Tint>(helper_quad_field_nearest(xI)));
}

template <typename T, typename Tring, int d>
inline void NearestInteger(QuadField<T, d> const &xI, QuadField<Tring, d> &xO) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  xO = QuadField<Tring, d>(
      UniversalScalarConversion<Tring, Tint>(helper_quad_field_nearest(xI)));
}

template <typename T, typename Tout, int d>
requires (!is_quad_field<Tout>::value)
inline void NearestInteger(QuadField<T, d> const &xI, Tout &xO) {
  using Tint = typename underlying_z_ring<T>::ring_type;
  xO = UniversalScalarConversion<Tout, Tint>(helper_quad_field_nearest(xI));
}

template <typename T, int d> struct is_totally_ordered<QuadField<T, d>> {
  static const bool value = true;
};

template <typename T, int d> struct is_ring_field<QuadField<T, d>> {
  static const bool value = is_ring_field<T>::value;
};

// FMA form (see is_fma_prefered). The compound form is fastest for
// QuadField: the QuadProd sinks accumulate through reused thread local
// scratches (allocation free after warm up), measured at 1.9x over the
// former fresh-temporary sinks and slightly ahead of the reused
// QuadField scratch.
template <typename T, int d> struct is_fma_prefered<QuadField<T, d>> {
  static const bool value = true;
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

// Fraction-free LU inverse is a clear win for quadratic fields over an exact
// base (benchmarked 2.5-6x vs classical); off over a floating-point base, and
// off over a ring base, where the inverse lies in the ring only when the
// determinant is a unit, so the generic non-field dispatch of Inverse -- go to
// the overlying field and come back -- is the correct behaviour.
template <typename T, int d>
struct use_fraction_free_lu<QuadField<T, d>> {
  static const bool value = is_exact_arithmetic<T>::value &&
                            is_ring_field<T>::value;
};

// Over a ring base there is no gcd of two ring elements to reduce a content
// with and no division to normalize with, so the vector canonicalization is
// done by ScalarCanonicalizationVectorRing below rather than through the
// overlying field.
template <typename T, int d>
struct has_ring_canonicalization<QuadField<T, d>> {
  static const bool value = !is_ring_field<T>::value;
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

// The other direction: a scalar that is not a quadratic field element enters
// the field as a + 0 sqrt(d). Always defined, no check needed, the base
// conversion alone can refuse. This is what lets a computation carry an
// integral quantity into a Q(sqrt(d)) matrix, the LLL size reduction of a
// form over the field being the first user: the reduction coefficient is a
// rational integer and has to be multiplied into the field.
template <typename T1, typename T2, int d>
requires (!is_quad_field<T1>::value)
inline void TYPE_CONVERSION(stc<T1> const &x1, QuadField<T2, d> &x2) {
  TYPE_CONVERSION(x1, x2.get_a());
  x2.get_b() = 0;
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

// The canonical representative of V on its ray, computed inside the ring.
// Z[sqrt(d)] has no gcd of two ring elements, so the content of V cannot be
// reduced away; what bounds the coefficients is the division by one entry,
// which the field does with an actual division. Here the same division is done
// through 1/s = conj(s) / N(s): V is multiplied by conj(s), which stays in the
// ring, and the content over the base ring of the result is divided out. The
// outcome is the primitive vector on the ray of V / s, which is what the field
// normalization produces too, but reached without a single fraction. The entry
// s is the one of smallest absolute value, as
// CanonicalizationSmallestCoefficientVectorPlusCoeff picks it.
template <typename T, int d>
MyVector<QuadField<T, d>>
ScalarCanonicalizationVectorRing(MyVector<QuadField<T, d>> const &V) {
  using Tquad = QuadField<T, d>;
  int n = V.size();
  int i_sma = -1;
  for (int i = 0; i < n; i++) {
    if (V(i) != 0) {
      if (i_sma == -1 || T_abs(V(i)) < T_abs(V(i_sma)))
        i_sma = i;
    }
  }
  if (i_sma == -1) {
    // The zero vector, already canonical.
    return V;
  }
  // The divisor is the absolute value of the smallest entry, so that the
  // direction of the result is the one the field normalization gives.
  Tquad s = T_abs(V(i_sma));
  Tquad conj_s(s.get_const_a(), -s.get_const_b());
  T norm = s.get_const_a() * s.get_const_a() -
           d * s.get_const_b() * s.get_const_b();
  // The norm carries the sign of the division by s, and the direction of the
  // result must be that of V / s.
  bool negate = (norm < 0);
  MyVector<Tquad> W(n);
  T g(0);
  for (int i = 0; i < n; i++) {
    Tquad val = V(i) * conj_s;
    if (negate)
      val = -val;
    g = GcdPair(g, val.get_const_a());
    g = GcdPair(g, val.get_const_b());
    W(i) = val;
  }
  if (g < 0)
    g = -g;
  if (g == 0 || g == 1)
    return W;
  MyVector<Tquad> Wred(n);
  for (int i = 0; i < n; i++)
    Wred(i) = Tquad(W(i).get_const_a() / g, W(i).get_const_b() / g);
  return Wred;
}

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYQUADFIELD_H_
// clang-format on
