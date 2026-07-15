// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_JET_NUMBER_H_
#define SRC_NUMBER_JET_NUMBER_H_

// clang-format off
#include "ExceptionsFunc.h"
#include "TemplateTraits.h"
#include "BasicNumberTypes.h"
#include <array>
#include <istream>
#include <ostream>
#include <utility>
// clang-format on

// A jet number: a truncated Taylor expansion c0 + c1 t + ... + cN t^N in an
// infinitesimal parameter t, kept at a fixed compile-time order N. This is a
// full numeric TYPE meant to be used as the scalar of generic templated code
// (matrices, determinants, CVP, ...): a matrix Q + t H over jet<T, N> propagates
// the order-N expansion of every derived quantity in a single exact
// computation, so the derivatives fall out of the coefficients (jet_deriv
// below) instead of being interpolated from many concrete-t samples.
//
// Every operation truncates back to order N, so a chain of operations stays at
// order N (truncation is a ring homomorphism onto T[t]/(t^{N+1})).
//
// ORDER (the total order making t an infinitesimal positive): the sign of a jet
// as t -> 0^+ is the sign of its first non-zero coefficient, so a is >= 0 iff
// its first non-zero coefficient is >= 0 (an all-zero jet is 0 >= 0). All the
// comparison operators derive from this leading-coefficient rule applied to the
// difference a - b. This is exactly the ordering needed for the geometric
// predicates (positivity, CVP inequalities, pivoting) on the segment t -> 0^+.

#ifdef SANITY_CHECK
#define SANITY_CHECK_JET_NUMBER
#endif

template <typename T, int N> class jet;

// Lazy product of two jets -- a minimal expression template (see the analogous
// RatProd / QuadProd). `a * b` returns this proxy (no computation, no
// allocation) instead of a fresh jet. The fast sinks evaluate the convolution
// directly into their own coefficient buffers, so no temporary jet (N+1 fresh T
// coefficients) is allocated per product:
//   prod  = a * b;   -> jet::operator=(jetProd)    (convolve into prod.c[])
//   acc  += a * b;   -> jet::operator+=(jetProd)   (fused into acc.c[])
//   acc  -= a * b;   -> jet::operator-=(jetProd)
//   jet r = a * b;   -> jet(jetProd)               (fresh, as before)
// Every other use materializes it into a jet through the operators after the
// class, so results are identical to the eager version. Holds references:
// consume within the same full-expression, do not bind with `auto`.
template <typename T, int N> struct jetProd {
  jet<T, N> const &x;
  jet<T, N> const &y;
};

template <typename T, int N> class jet {
public:
  std::array<T, N + 1> c; // c[k] = coefficient of t^k, 0 <= k <= N

  // Zero.
  jet() { c.fill(T(0)); }
  // Constant jets (non-explicit so generic code's T(0), T(1), 0, ... promote).
  jet(int v) {
    c.fill(T(0));
    c[0] = T(v);
  }
  jet(T const &v) {
    c.fill(T(0));
    c[0] = v;
  }
  // Construct from a lazy product a*b: materialize the convolution into this
  // fresh jet. Also the implicit jetProd -> jet conversion for non-sink uses.
  jet(jetProd<T, N> const &e) {
    c.fill(T(0));
    T prod;
    for (int i = 0; i <= N; i++)
      for (int j = 0; i + j <= N; j++) {
        prod = e.x.c[i] * e.y.c[j];
        c[i + j] += prod;
      }
  }
  // The infinitesimal t itself (c1 = 1, all others 0).
  static jet var() {
    jet r;
    r.c[1] = T(1);
    return r;
  }
  // The jet whose coefficient array is exactly `coeffs`.
  static jet from_coeffs(std::array<T, N + 1> coeffs) {
    jet r;
    r.c = std::move(coeffs);
    return r;
  }

  T const &operator[](int k) const { return c[k]; }
  T &operator[](int k) { return c[k]; }

  // ---- additive group (member unary and compound assignment) ----
  jet operator-() const {
    jet r;
    for (int k = 0; k <= N; k++)
      r.c[k] = -c[k];
    return r;
  }
  jet &operator+=(jet const &o) {
    for (int k = 0; k <= N; k++)
      c[k] += o.c[k];
    return *this;
  }
  jet &operator-=(jet const &o) {
    for (int k = 0; k <= N; k++)
      c[k] -= o.c[k];
    return *this;
  }
  // Assign from a lazy product a*b: convolve directly into this->c[], reusing
  // its buffers (no fresh jet is allocated). If this aliases an operand, fall
  // back to a materialized temporary, since the convolution reads the operands'
  // coefficients while it would be overwriting this->c[].
  jet &operator=(jetProd<T, N> const &e) {
    if (this == &e.x || this == &e.y) {
      jet tmp(e);
      c = std::move(tmp.c);
      return *this;
    }
    for (int k = 0; k <= N; k++)
      c[k] = T(0);
    T prod;
    for (int i = 0; i <= N; i++)
      for (int j = 0; i + j <= N; j++) {
        prod = e.x.c[i] * e.y.c[j];
        c[i + j] += prod;
      }
    return *this;
  }
  // Fused accumulate of a lazy product: this += a*b, convolved directly into
  // this->c[] with no temporary jet (aliasing handled as in operator=).
  jet &operator+=(jetProd<T, N> const &e) {
    if (this == &e.x || this == &e.y) {
      jet tmp(e);
      for (int k = 0; k <= N; k++)
        c[k] += tmp.c[k];
      return *this;
    }
    T prod;
    for (int i = 0; i <= N; i++)
      for (int j = 0; i + j <= N; j++) {
        prod = e.x.c[i] * e.y.c[j];
        c[i + j] += prod;
      }
    return *this;
  }
  // Fused subtract of a lazy product: this -= a*b.
  jet &operator-=(jetProd<T, N> const &e) {
    if (this == &e.x || this == &e.y) {
      jet tmp(e);
      for (int k = 0; k <= N; k++)
        c[k] -= tmp.c[k];
      return *this;
    }
    T prod;
    for (int i = 0; i <= N; i++)
      for (int j = 0; i + j <= N; j++) {
        prod = e.x.c[i] * e.y.c[j];
        c[i + j] -= prod;
      }
    return *this;
  }
  jet &operator*=(jet const &o) {
    *this = (*this) * o;
    return *this;
  }
  jet &operator/=(jet const &o) {
    *this = (*this) / o;
    return *this;
  }

  // Index of the first non-zero coefficient, or N+1 if the jet is zero.
  int leading_index() const {
    for (int k = 0; k <= N; k++)
      if (c[k] != T(0))
        return k;
    return N + 1;
  }
  // Sign as t -> 0^+ : sign of the first non-zero coefficient (0 if all zero).
  int sign() const {
    int k = leading_index();
    if (k > N)
      return 0;
    if (c[k] > T(0))
      return 1;
    return -1;
  }

  // Evaluate the truncated series at a concrete value t (Horner). Used by the
  // higher-level cross-validation code.
  T eval(T const &t) const {
    T s = c[N];
    for (int k = N - 1; k >= 0; k--)
      s = s * t + c[k];
    return s;
  }

  // ---- ring operations and total order, as hidden friends so that a mixed
  // expression like `x != 0` or `alpha * x` (with an int / T literal) promotes
  // the literal to a constant jet. Plain function templates would not apply the
  // implicit conversion during argument deduction. ----
  friend jet operator+(jet a, jet const &b) {
    a += b;
    return a;
  }
  friend jet operator-(jet a, jet const &b) {
    a -= b;
    return a;
  }
  // Lazy: returns a jetProd proxy (see above), evaluated in place by the
  // consumer (operator= / operator+= / operator-= / the constructor).
  friend jetProd<T, N> operator*(jet const &a, jet const &b) {
    return jetProd<T, N>{a, b};
  }

  // Inverse 1/f of a jet with non-zero constant term, from f * f^{-1} = 1:
  // b_0 = 1/c_0 and b_k = -b_0 sum_{j=1}^{k} c_j b_{k-j}.
  friend jet inverse(jet const &f) {
#ifdef SANITY_CHECK_JET_NUMBER
    if (f.c[0] == T(0)) {
      std::cerr << "jet_number: inverse of a jet with zero constant term "
                   "(would be a Laurent series, unsupported)\n";
      throw TerminalException{1};
    }
#endif
    jet r;
    r.c[0] = T(1) / f.c[0];
    // Hoisted scratch (buffer reuse across iterations, as in operator*).
    T s, prod;
    for (int k = 1; k <= N; k++) {
      s = T(0);
      for (int j = 1; j <= k; j++) {
        prod = f.c[j] * r.c[k - j];
        s += prod;
      }
      prod = r.c[0] * s;
      r.c[k] = -prod;
    }
#ifdef SANITY_CHECK_JET_NUMBER
    // Cheap invariant: f * f^{-1} == 1 (one order-N convolution).
    jet chk = f * r;
    if (chk.c[0] != T(1)) {
      std::cerr << "jet_number: inverse check failed (c0 != 1)\n";
      throw TerminalException{1};
    }
    for (int k = 1; k <= N; k++)
      if (chk.c[k] != T(0)) {
        std::cerr << "jet_number: inverse check failed (c" << k << " != 0)\n";
        throw TerminalException{1};
      }
#endif
    return r;
  }

  // Division a / b, defined only when b is a unit (non-zero constant term), in
  // which case it is a * inverse(b). Division by a jet with c0 == 0 (a
  // zero-divisor of T[t]/(t^{N+1})) has no representable quotient in general
  // (the true result is a Laurent series). Such a division must never occur: the
  // computations that would otherwise divide by a rationally-degenerate quantity
  // (e.g. a simplex volume alpha*t^d that vanishes at t = 0) are routed through
  // the division-free determinant (see determinant_division_free below) instead.
  // Hitting this with a non-unit b is therefore a programming error, flagged
  // under SANITY_CHECK rather than papered over.
  friend jet operator/(jet const &a, jet const &b) {
#ifdef SANITY_CHECK_JET_NUMBER
    if (b.c[0] == T(0)) {
      std::cerr << "jet_number: division by a jet with zero constant term "
                   "(zero-divisor); the caller should use a division-free path\n";
      throw TerminalException{1};
    }
#endif
    return a * inverse(b);
  }

  // Total order by the leading coefficient of the difference.
  friend bool operator==(jet const &a, jet const &b) {
    for (int k = 0; k <= N; k++)
      if (a.c[k] != b.c[k])
        return false;
    return true;
  }
  friend bool operator!=(jet const &a, jet const &b) { return !(a == b); }
  friend bool operator<(jet const &a, jet const &b) {
    return (a - b).sign() < 0;
  }
  friend bool operator>(jet const &a, jet const &b) {
    return (a - b).sign() > 0;
  }
  friend bool operator<=(jet const &a, jet const &b) {
    return (a - b).sign() <= 0;
  }
  friend bool operator>=(jet const &a, jet const &b) {
    return (a - b).sign() >= 0;
  }
};

// ---------------------------------------------------------------------------
// jetProd (the lazy a*b proxy) as a first-class value. Every use other than the
// in-place sinks above materializes the proxy into a jet and delegates to the
// ordinary jet operators, so results are identical to the eager version. The
// arithmetic operators return jet explicitly so that a jetProd produced on the
// right-hand side is materialized before its operand temporaries die.
// ---------------------------------------------------------------------------
template <typename T, int N>
inline jet<T, N> const &jet_eval(jet<T, N> const &x) {
  return x;
}
template <typename T, int N>
inline jet<T, N> jet_eval(jetProd<T, N> const &e) {
  return jet<T, N>(e);
}

#define JET_JETPROD_ARITH(OP)                                                  \
  template <typename T, int N>                                                 \
  inline jet<T, N> operator OP(jetProd<T, N> const &a,                         \
                               jetProd<T, N> const &b) {                       \
    return jet_eval(a) OP jet_eval(b);                                         \
  }                                                                            \
  template <typename T, int N>                                                 \
  inline jet<T, N> operator OP(jetProd<T, N> const &a, jet<T, N> const &b) {   \
    return jet_eval(a) OP b;                                                   \
  }                                                                            \
  template <typename T, int N>                                                 \
  inline jet<T, N> operator OP(jet<T, N> const &a, jetProd<T, N> const &b) {   \
    return a OP jet_eval(b);                                                   \
  }
JET_JETPROD_ARITH(+)
JET_JETPROD_ARITH(-)
JET_JETPROD_ARITH(*)
JET_JETPROD_ARITH(/)
#undef JET_JETPROD_ARITH

#define JET_JETPROD_CMP(OP)                                                    \
  template <typename T, int N>                                                 \
  inline bool operator OP(jetProd<T, N> const &a, jetProd<T, N> const &b) {    \
    return jet_eval(a) OP jet_eval(b);                                         \
  }                                                                            \
  template <typename T, int N>                                                 \
  inline bool operator OP(jetProd<T, N> const &a, jet<T, N> const &b) {        \
    return jet_eval(a) OP b;                                                   \
  }                                                                            \
  template <typename T, int N>                                                 \
  inline bool operator OP(jet<T, N> const &a, jetProd<T, N> const &b) {        \
    return a OP jet_eval(b);                                                   \
  }
JET_JETPROD_CMP(==)
JET_JETPROD_CMP(!=)
JET_JETPROD_CMP(<)
JET_JETPROD_CMP(>)
JET_JETPROD_CMP(<=)
JET_JETPROD_CMP(>=)
#undef JET_JETPROD_CMP

template <typename T, int N>
inline jet<T, N> operator-(jetProd<T, N> const &e) {
  return -jet_eval(e);
}

// The constant term (value at t = 0) of a jet is c0. This overloads the generic
// constant_term (the identity, in TemplateTraits.h) so that the combinatorial /
// canonical-form subroutines of a scalar-templated computation (which need a
// concrete field element) recover the t = 0 data from a jet Gram matrix, while
// the numeric parts keep the full expansion.
template <typename T, int N> T const &constant_term(jet<T, N> const &j) {
  return j.c[0];
}

// The k-th derivative at t = 0: k! times the coefficient of t^k.
template <typename T, int N> T jet_deriv(jet<T, N> const &j, int k) {
  T fact(1);
  for (int i = 2; i <= k; i++)
    fact *= T(i);
  return fact * j.c[k];
}

// ---- number-type traits so jet<T, N> flows through the generic matrix code
// (MyMatrix, DeterminantMat, Inverse, SolutionMat, ...). The truncated series
// ring is a local ring, not a field, but the hand-written Gaussian-elimination
// kernels behave as if over a field as long as the pivots have non-zero
// constant term -- which the pivot-cost rule below guarantees.
template <typename T, int N> struct is_ring_field<jet<T, N>> {
  static const bool value = is_ring_field<T>::value;
};
// FMA form (see is_fma_prefered). The direct/fused form is fastest for jet
// (measured): operator+=(jetProd) / operator-=(jetProd) convolve the product
// directly into the accumulator in a single pass, with no temporary jet.
template <typename T, int N> struct is_fma_prefered<jet<T, N>> {
  static const bool value = true;
};
template <typename T, int N> struct is_totally_ordered<jet<T, N>> {
  static const bool value = true;
};
template <typename T, int N> struct is_exact_arithmetic<jet<T, N>> {
  static const bool value = is_exact_arithmetic<T>::value;
};
template <typename T, int N> struct overlying_field<jet<T, N>> {
  typedef jet<typename overlying_field<T>::field_type, N> field_type;
};
template <typename T, int N> struct underlying_ring<jet<T, N>> {
  typedef jet<typename underlying_ring<T>::ring_type, N> ring_type;
};
template <typename T, int N> struct is_implementation_of_Q<jet<T, N>> {
  static const bool value = false;
};
template <typename T, int N> struct is_implementation_of_Z<jet<T, N>> {
  static const bool value = false;
};
template <typename T, int N> struct underlying_totally_ordered_ring<jet<T, N>> {
  typedef jet<typename underlying_totally_ordered_ring<T>::real_type, N>
      real_type;
};
// The truncated jet ring has zero divisors (any jet with c0 == 0), so an
// elimination-based determinant can be forced to divide by one. Route jet
// determinants through the division-free dispatch (see DeterminantMat).
template <typename T, int N> struct determinant_division_free<jet<T, N>> {
  static const bool value = true;
};

// A jet represents an integer iff it is a constant with integer constant term.
template <typename T, int N> bool IsInteger(jet<T, N> const &x) {
  for (int k = 1; k <= N; k++)
    if (x.c[k] != T(0))
      return false;
  return IsInteger(x.c[0]);
}

// Conversion to double goes through the constant term (the value at t = 0); it
// is only used for the floating-point diagnostics.
template <typename T, int N>
inline void TYPE_CONVERSION(stc<jet<T, N>> const &x1, double &x2) {
  stc<T> a1{x1.val.c[0]};
  TYPE_CONVERSION(a1, x2);
}

// Conversion of a scalar (a rational / integer at t = 0) into a jet: the value
// becomes the constant term, all higher coefficients zero. This is how the
// t = 0 data (SHV, EXT, the Gram matrix Q) is lifted into the jet computation.
template <typename Tin, typename T, int N>
inline void TYPE_CONVERSION(stc<Tin> const &x1, jet<T, N> &x2) {
  T val;
  TYPE_CONVERSION(x1, val);
  x2 = jet<T, N>(val);
}
// Jet -> jet of the same type is the identity (more specialized than the scalar
// lift above, so it wins for UniversalScalarConversion<jet, jet>).
template <typename T, int N>
inline void TYPE_CONVERSION(stc<jet<T, N>> const &x1, jet<T, N> &x2) {
  x2 = x1.val;
}

// Reading a jet from a stream (used by ParseScalar) fills the constant term.
template <typename T, int N>
std::istream &operator>>(std::istream &is, jet<T, N> &j) {
  j = jet<T, N>();
  is >> j.c[0];
  return is;
}

// Pivot cost of a jet (see PivotCost.h): the order (index of the first non-zero
// coefficient) together with the base-field cost of that leading coefficient.
// Smallest order is preferred, so the elimination always selects a pivot with
// non-zero constant term (order 0) when one exists -- exactly the invertible
// pivots, which for Q + t H with Q non-singular is guaranteed at every step.
// Ties on the order are broken by the base-field rule on the leading
// coefficient.
//
// The cost is a dedicated struct in the global namespace (not std::pair) so that
// the generic SelectBestPivot in the matrix kernels finds is_preferable_pivot by
// argument-dependent lookup even though that overload lives here, after the
// kernels are defined.
template <typename C> struct JetPivotCost {
  size_t order;
  C base;
};
template <typename T, int N>
JetPivotCost<decltype(f_cost_pivot(std::declval<T>()))>
f_cost_pivot(jet<T, N> const &j) {
  int k = j.leading_index();
  if (k > N)
    k = N; // an all-zero jet is never a pivot; guard the coefficient access
  return {static_cast<size_t>(k), f_cost_pivot(j.c[k])};
}
template <typename C>
bool is_preferable_pivot(JetPivotCost<C> const &x, JetPivotCost<C> const &y) {
  if (x.order != y.order)
    return x.order < y.order;
  return is_preferable_pivot(x.base, y.base);
}

template <typename T, int N>
std::ostream &operator<<(std::ostream &os, jet<T, N> const &j) {
  os << "[";
  for (int k = 0; k <= N; k++) {
    if (k > 0)
      os << ", ";
    os << j.c[k];
  }
  os << "]";
  return os;
}

// Hash over all coefficients (so distinct jets -- including two that share a
// constant term but differ at higher order -- hash apart, as needed by the
// weight-matrix / canonical-form coloring that resolves the t -> 0^+ structure).
namespace std {
template <typename T, int N> struct hash<jet<T, N>> {
  std::size_t operator()(const jet<T, N> &x) const {
    auto combine_hash = [](size_t &seed, size_t new_hash) -> void {
      seed ^= new_hash + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    };
    size_t seed = std::hash<int>()(N);
    for (int k = 0; k <= N; k++)
      combine_hash(seed, std::hash<T>()(x.c[k]));
    return seed;
  }
};
// clang-format off
}  // namespace std
// clang-format on

// clang-format off
#endif  // SRC_NUMBER_JET_NUMBER_H_
// clang-format on
