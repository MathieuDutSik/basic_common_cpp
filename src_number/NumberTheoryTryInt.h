// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYTRYINT_H_
#define SRC_NUMBER_NUMBERTHEORYTRYINT_H_

// clang-format off
#include "TypeConversion.h"
#include "rational.h"
#include <cstdint>
#include <iostream>
#include <limits>
#include <type_traits>
// clang-format on

//
// TryInt64: an int64_t with DEFERRED overflow detection.
//
// SafeInt64 throws at the very operation that overflows, which puts a branch
// to a throwing path behind every single +, -, *. That inhibits instruction
// level parallelism and SIMD vectorization in tight loops (measured on the
// Bareiss determinant kernel: SafeInt64 is ~12x over raw int64_t, TryInt64
// ~1.6x). TryInt64 instead runs the plain machine operation, records an
// overflow in the thread-local flag is_correct, and lets the computation
// continue with wrapped values. Strategically placed calls to
// terminate_in_arithmetic_error<T>() -- out of the tight loop but not
// completely out of scope, typically once per outer iteration -- check the
// flag and throw TryIntException, so a wrapped value never escapes into a
// result. The caller converts the input from an implementation-of-Z type T
// (throwing TryIntException when an entry does not fit), runs the expensive
// code over TryInt64, and on TryIntException falls back to running it over
// the original T, which is slower but always correct.
//
// The divisions guard against a zero or overflowing divisor (a divisor can be
// a wrapped garbage zero after an earlier overflow) by flagging is_correct
// instead of trapping, so the deferred checking never crashes in between two
// terminate_in_arithmetic_error calls.
//

// The inline (single entity program-wide) is required: the inline operators
// below store into it, so a per-translation-unit static would be an ODR trap.
inline thread_local bool is_correct = true;

inline int64_t try_add_int64(int64_t a, int64_t b) {
  int64_t r;
  is_correct &= !__builtin_add_overflow(a, b, &r);
  return r;
}

inline int64_t try_subtract_int64(int64_t a, int64_t b) {
  int64_t r;
  is_correct &= !__builtin_sub_overflow(a, b, &r);
  return r;
}

inline int64_t try_multiply_int64(int64_t a, int64_t b) {
  int64_t r;
  is_correct &= !__builtin_mul_overflow(a, b, &r);
  return r;
}

inline int64_t try_negate_int64(int64_t a) { return try_subtract_int64(0, a); }

inline int64_t try_divide_int64(int64_t a, int64_t b) {
  if (b == 0 || (b == -1 && a == std::numeric_limits<int64_t>::min())) {
    is_correct = false;
    return 0;
  }
  return a / b;
}

inline int64_t try_modulo_int64(int64_t a, int64_t b) {
  if (b == 0 || (b == -1 && a == std::numeric_limits<int64_t>::min())) {
    is_correct = false;
    return 0;
  }
  return a % b;
}

struct TryInt64 {
public:
  using Tint = int64_t;

private:
  Tint val;

public:
  TryInt64() : val(0) {}
  // Single constructor accepting any built-in integral type; see the
  // corresponding SafeInt64 constructor for the rationale.
  template <class U,
            std::enable_if_t<std::is_integral_v<U>, int> = 0>
  TryInt64(U const &x) : val(static_cast<Tint>(x)) {}
  TryInt64(TryInt64 const &x) : val(x.val) {}
  double get_d() const { return static_cast<double>(val); }
  TryInt64 &operator=(TryInt64 const &u) {
    val = u.val;
    return *this;
  }
  Tint &get_val() { return val; }
  const Tint &get_const_val() const { return val; }
  TryInt64 &operator+=(TryInt64 const &x) {
    val = try_add_int64(val, x.val);
    return *this;
  }
  TryInt64 &operator-=(TryInt64 const &x) {
    val = try_subtract_int64(val, x.val);
    return *this;
  }
  friend TryInt64 operator+(TryInt64 const &x, TryInt64 const &y) {
    TryInt64 z;
    z.val = try_add_int64(x.val, y.val);
    return z;
  }
  friend TryInt64 operator+(Tint const &x, TryInt64 const &y) {
    TryInt64 z;
    z.val = try_add_int64(x, y.val);
    return z;
  }
  friend TryInt64 operator+(TryInt64 const &x, Tint const &y) {
    TryInt64 z;
    z.val = try_add_int64(x.val, y);
    return z;
  }
  friend TryInt64 operator-(TryInt64 const &x, TryInt64 const &y) {
    TryInt64 z;
    z.val = try_subtract_int64(x.val, y.val);
    return z;
  }
  friend TryInt64 operator-(TryInt64 const &x) {
    TryInt64 z;
    z.val = try_negate_int64(x.val);
    return z;
  }
  TryInt64 &operator++() {
    val = try_add_int64(val, 1);
    return *this;
  }
  TryInt64 operator++(int) {
    TryInt64 tmp = *this;
    val = try_add_int64(val, 1);
    return tmp;
  }
  TryInt64 &operator--() {
    val = try_subtract_int64(val, 1);
    return *this;
  }
  TryInt64 operator--(int) {
    TryInt64 tmp = *this;
    val = try_subtract_int64(val, 1);
    return tmp;
  }
  TryInt64 &operator*=(TryInt64 const &x) {
    val = try_multiply_int64(val, x.val);
    return *this;
  }
  friend TryInt64 operator*(TryInt64 const &x, TryInt64 const &y) {
    TryInt64 z;
    z.val = try_multiply_int64(x.val, y.val);
    return z;
  }
  friend TryInt64 operator/(TryInt64 const &x, TryInt64 const &y) {
    TryInt64 z;
    z.val = try_divide_int64(x.val, y.val);
    return z;
  }
  friend TryInt64 operator%(TryInt64 const &x, TryInt64 const &y) {
    TryInt64 z;
    z.val = try_modulo_int64(x.val, y.val);
    return z;
  }
  friend TryInt64 operator*(Tint const &x, TryInt64 const &y) {
    TryInt64 z;
    z.val = try_multiply_int64(x, y.val);
    return z;
  }
  friend TryInt64 operator*(TryInt64 const &x, Tint const &y) {
    TryInt64 z;
    z.val = try_multiply_int64(x.val, y);
    return z;
  }
  friend std::ostream &operator<<(std::ostream &os, TryInt64 const &v) {
    return os << v.val;
  }
  friend std::istream &operator>>(std::istream &is, TryInt64 &v) {
    is >> v.val;
    return is;
  }
  //
  friend bool operator==(TryInt64 const &x, TryInt64 const &y) {
    return x.val == y.val;
  }
  friend bool operator==(TryInt64 const &x, int const &y) {
    return x.val == static_cast<Tint>(y);
  }
  //
  friend bool operator!=(TryInt64 const &x, TryInt64 const &y) {
    return x.val != y.val;
  }
  friend bool operator!=(TryInt64 const &x, int const &y) {
    return x.val != static_cast<Tint>(y);
  }
  //
  friend bool operator>=(TryInt64 const &x, TryInt64 const &y) {
    return x.val >= y.val;
  }
  friend bool operator>=(TryInt64 const &x, int const &y) {
    return x.val >= static_cast<Tint>(y);
  }
  //
  friend bool operator<=(TryInt64 const &x, TryInt64 const &y) {
    return x.val <= y.val;
  }
  friend bool operator<=(TryInt64 const &x, int const &y) {
    return x.val <= static_cast<Tint>(y);
  }
  //
  friend bool operator>(TryInt64 const &x, TryInt64 const &y) {
    return x.val > y.val;
  }
  friend bool operator>(TryInt64 const &x, int const &y) {
    return x.val > static_cast<Tint>(y);
  }
  //
  friend bool operator<(TryInt64 const &x, TryInt64 const &y) {
    return x.val < y.val;
  }
  friend bool operator<(TryInt64 const &x, int const &y) {
    return x.val < static_cast<Tint>(y);
  }
};

// The deferred overflow check. Trivial for every type; the deferred-checking
// types check the is_correct flag, reset it, and throw TryIntException so the
// caller can fall back to an exact type. Call it out of the tight loop but not
// completely out of scope (typically once per outer loop iteration) so that a
// computation on wrapped values ends prematurely.
template <typename T> inline void terminate_in_arithmetic_error() {}

template <> inline void terminate_in_arithmetic_error<TryInt64>() {
  if (!is_correct) {
    is_correct = true;
    throw TryIntException{1};
  }
}

template <> inline void terminate_in_arithmetic_error<Rational<TryInt64>>() {
  terminate_in_arithmetic_error<TryInt64>();
}

// Conversion from an implementation-of-Z type T (e.g. mpz_class) into
// TryInt64, throwing TryIntException when the value does not fit into
// int64_t. Some TYPE_CONVERSION into int64_t truncate silently (e.g.
// mpz_class via get_si), so the conversion is verified by a round trip. The
// scratch holds the round-trip value: passing it in from the caller keeps
// its allocation (e.g. the mpz limb) out of entrywise conversion loops.
// Number types with a cheap exact fit test provide a non-template overload
// next to their own definition (e.g. mpz_class in NumberTheoryGmp.h), which
// overload resolution prefers to this generic version.
template <typename T>
TryInt64 ConvertToTryInt64(T const &val, T &scratch) {
  int64_t v;
  try {
    TYPE_CONVERSION(stc<T>{val}, v);
  } catch (ConversionException const &) {
    throw TryIntException{1};
  }
  TYPE_CONVERSION(stc<int64_t>{v}, scratch);
  if (!(scratch == val)) {
    throw TryIntException{1};
  }
  return TryInt64(v);
}

template <typename T> TryInt64 ConvertToTryInt64(T const &val) {
  T scratch;
  return ConvertToTryInt64(val, scratch);
}

template <typename T> T ConvertFromTryInt64(TryInt64 const &val) {
  int64_t const &v = val.get_const_val();
  T ret;
  TYPE_CONVERSION(stc<int64_t>{v}, ret);
  return ret;
}

inline void TYPE_CONVERSION(stc<int64_t> const &a1, TryInt64 &a2) {
  a2 = TryInt64(a1.val);
}

inline void TYPE_CONVERSION(stc<TryInt64> const &a1, int64_t &a2) {
  a2 = a1.val.get_const_val();
}

// Traits: same profile as SafeInt64.

template <> struct is_implementation_of_Z<TryInt64> {
  static const bool value = true;
};

template <> struct is_implementation_of_Q<TryInt64> {
  static const bool value = false;
};

template <> struct is_euclidean_domain<TryInt64> {
  static const bool value = true;
};

template <> struct is_ring_field<TryInt64> {
  static const bool value = false;
};

template <> struct is_totally_ordered<TryInt64> {
  static const bool value = true;
};

// Exact provided the terminate_in_arithmetic_error discipline is followed: a
// wrapped value is always reported before a result is returned.
template <> struct is_exact_arithmetic<TryInt64> {
  static const bool value = true;
};

// Native-width integer with no heap allocation: fused form, nothing for a
// scratch to save (see is_fma_prefered).
template <> struct is_fma_prefered<TryInt64> {
  static const bool value = true;
};

template <> struct underlying_ring<TryInt64> {
  typedef TryInt64 ring_type;
};

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYTRYINT_H_
// clang-format on
