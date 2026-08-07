// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_NUMBERTHEORYTRYINT_H_
#define SRC_NUMBER_NUMBERTHEORYTRYINT_H_

// clang-format off
#include "QuoIntFcts.h"
#include "ResidueQuotient.h"
#include "TypeConversion.h"
#include "rational.h"
#include <cstdint>
#include <iostream>
#include <limits>
#include <type_traits>
// clang-format on

//
// Deferred-overflow int64 arithmetic, in TWO flavours sharing one design.
//
// SafeInt64 throws at the very operation that overflows, which puts a branch
// to a throwing path behind every single +, -, *. That inhibits instruction
// level parallelism and SIMD vectorization in tight loops. The types here
// instead run the plain machine operation, record an overflow in the
// thread-local flag is_correct, and let the computation continue with wrapped
// values. Strategically placed calls to terminate_in_arithmetic_error<T>() --
// out of the tight loop but not completely out of scope, typically once per
// outer iteration -- check the flag and throw TryIntException, so a wrapped
// value never escapes into a result. The caller converts the input from an
// implementation-of-Z type T (throwing TryIntException when an entry does not
// fit), runs the expensive code over the try-type, and on TryIntException
// falls back to running it over the original T, which is slower but correct.
//
// The two flavours differ ONLY in how an overflow is detected:
//
//  --- TryCarryInt64 (policy try_carry_policy): exact detection via the
//      hardware carry/overflow flag (__builtin_*_overflow). It never
//      false-flags, so it falls back only on a genuine overflow, but each
//      op does a read-modify-write of is_correct (`is_correct &= ...`),
//      a per-op serializing store that blocks vectorization. Best when the
//      arithmetic is not the bottleneck or coefficients are medium-sized.
//
//  --- TrySimdInt64 (policy try_simd_policy): a-priori PREDICTION, the trick
//      lrslib's SAFE arithmetic uses. It compares the OPERANDS against a
//      conservative safe threshold (a product fits when |a|,|b| <=
//      floor(sqrt(2^63-1)); a sum/difference fits when |a|,|b| <= 2^62-1);
//      inside the zone the exact result provably fits. Two wins: the check
//      is a plain compare, not a widening multiply + flag read; and on the
//      common (safe) path NOTHING is written to is_correct, so there is no
//      per-op serializing store and the loop can vectorize (SIMD friendly).
//      The price is that the prediction is conservative: a large operand
//      times a tiny one may be flagged even though the exact result fits,
//      causing a more frequent fallback to the exact type.
//
// Both flavours share the thread-local is_correct flag (a computation uses
// one flavour at a time), the division guards, the terminate check, the
// conversions and all the traits. TryInt64 is an alias for TryCarryInt64,
// the exact-detection default, so existing call sites are unchanged.
//
// The divisions guard against a zero or overflowing divisor (a divisor can be
// a wrapped garbage value after an earlier overflow) by flagging is_correct
// instead of trapping, so the deferred checking never crashes in between two
// terminate_in_arithmetic_error calls. The out-of-zone arithmetic in
// try_simd_policy is done in uint64_t so a flagged operation wraps with
// defined behavior rather than signed-overflow UB; the wrapped value is
// garbage but is never consumed.
//

// The inline (single entity program-wide) is required: the inline operators
// below store into it, so a per-translation-unit static would be an ODR trap.
inline thread_local bool is_correct = true;

// Exact carry-flag detection.
struct try_carry_policy {
  static inline int64_t add(int64_t a, int64_t b) {
    int64_t r;
    is_correct &= !__builtin_add_overflow(a, b, &r);
    return r;
  }
  static inline int64_t sub(int64_t a, int64_t b) {
    int64_t r;
    is_correct &= !__builtin_sub_overflow(a, b, &r);
    return r;
  }
  static inline int64_t mul(int64_t a, int64_t b) {
    int64_t r;
    is_correct &= !__builtin_mul_overflow(a, b, &r);
    return r;
  }
};

// Conservative a-priori prediction, no carry-flag dependency, SIMD friendly.
struct try_simd_policy {
  static constexpr int64_t mul_safe = 3037000499LL;        // floor(sqrt(2^63-1))
  static constexpr int64_t add_safe = 4611686018427387903LL; // 2^62 - 1
  static inline int64_t add(int64_t a, int64_t b) {
    if (a > add_safe || a < -add_safe || b > add_safe || b < -add_safe)
      is_correct = false;
    return static_cast<int64_t>(static_cast<uint64_t>(a) +
                                static_cast<uint64_t>(b));
  }
  static inline int64_t sub(int64_t a, int64_t b) {
    if (a > add_safe || a < -add_safe || b > add_safe || b < -add_safe)
      is_correct = false;
    return static_cast<int64_t>(static_cast<uint64_t>(a) -
                                static_cast<uint64_t>(b));
  }
  static inline int64_t mul(int64_t a, int64_t b) {
    if (a > mul_safe || a < -mul_safe || b > mul_safe || b < -mul_safe)
      is_correct = false;
    return static_cast<int64_t>(static_cast<uint64_t>(a) *
                                static_cast<uint64_t>(b));
  }
};

// Division and modulo are policy-independent: the guard is the same in both
// flavours (a garbage divisor after an earlier overflow must flag, not trap).
inline int64_t tryint_divide_int64(int64_t a, int64_t b) {
  if (b == 0 || (b == -1 && a == std::numeric_limits<int64_t>::min())) {
    is_correct = false;
    return 0;
  }
  return a / b;
}

inline int64_t tryint_modulo_int64(int64_t a, int64_t b) {
  if (b == 0 || (b == -1 && a == std::numeric_limits<int64_t>::min())) {
    is_correct = false;
    return 0;
  }
  return a % b;
}

template <class Pol> struct TryIntGen {
public:
  using Tint = int64_t;

private:
  Tint val;

public:
  TryIntGen() : val(0) {}
  // Single constructor accepting any built-in integral type; see the
  // corresponding SafeInt64 constructor for the rationale.
  template <class U, std::enable_if_t<std::is_integral_v<U>, int> = 0>
  TryIntGen(U const &x) : val(static_cast<Tint>(x)) {}
  TryIntGen(TryIntGen const &x) : val(x.val) {}
  double get_d() const { return static_cast<double>(val); }
  TryIntGen &operator=(TryIntGen const &u) {
    val = u.val;
    return *this;
  }
  Tint &get_val() { return val; }
  const Tint &get_const_val() const { return val; }
  TryIntGen &operator+=(TryIntGen const &x) {
    val = Pol::add(val, x.val);
    return *this;
  }
  TryIntGen &operator-=(TryIntGen const &x) {
    val = Pol::sub(val, x.val);
    return *this;
  }
  friend TryIntGen operator+(TryIntGen const &x, TryIntGen const &y) {
    TryIntGen z;
    z.val = Pol::add(x.val, y.val);
    return z;
  }
  friend TryIntGen operator+(Tint const &x, TryIntGen const &y) {
    TryIntGen z;
    z.val = Pol::add(x, y.val);
    return z;
  }
  friend TryIntGen operator+(TryIntGen const &x, Tint const &y) {
    TryIntGen z;
    z.val = Pol::add(x.val, y);
    return z;
  }
  friend TryIntGen operator-(TryIntGen const &x, TryIntGen const &y) {
    TryIntGen z;
    z.val = Pol::sub(x.val, y.val);
    return z;
  }
  friend TryIntGen operator-(TryIntGen const &x) {
    TryIntGen z;
    z.val = Pol::sub(0, x.val);
    return z;
  }
  TryIntGen &operator++() {
    val = Pol::add(val, 1);
    return *this;
  }
  TryIntGen operator++(int) {
    TryIntGen tmp = *this;
    val = Pol::add(val, 1);
    return tmp;
  }
  TryIntGen &operator--() {
    val = Pol::sub(val, 1);
    return *this;
  }
  TryIntGen operator--(int) {
    TryIntGen tmp = *this;
    val = Pol::sub(val, 1);
    return tmp;
  }
  TryIntGen &operator*=(TryIntGen const &x) {
    val = Pol::mul(val, x.val);
    return *this;
  }
  friend TryIntGen operator*(TryIntGen const &x, TryIntGen const &y) {
    TryIntGen z;
    z.val = Pol::mul(x.val, y.val);
    return z;
  }
  friend TryIntGen operator/(TryIntGen const &x, TryIntGen const &y) {
    TryIntGen z;
    z.val = tryint_divide_int64(x.val, y.val);
    return z;
  }
  friend TryIntGen operator%(TryIntGen const &x, TryIntGen const &y) {
    TryIntGen z;
    z.val = tryint_modulo_int64(x.val, y.val);
    return z;
  }
  friend TryIntGen operator*(Tint const &x, TryIntGen const &y) {
    TryIntGen z;
    z.val = Pol::mul(x, y.val);
    return z;
  }
  friend TryIntGen operator*(TryIntGen const &x, Tint const &y) {
    TryIntGen z;
    z.val = Pol::mul(x.val, y);
    return z;
  }
  friend std::ostream &operator<<(std::ostream &os, TryIntGen const &v) {
    return os << v.val;
  }
  friend std::istream &operator>>(std::istream &is, TryIntGen &v) {
    is >> v.val;
    return is;
  }
  //
  friend bool operator==(TryIntGen const &x, TryIntGen const &y) {
    return x.val == y.val;
  }
  friend bool operator==(TryIntGen const &x, int const &y) {
    return x.val == static_cast<Tint>(y);
  }
  //
  friend bool operator!=(TryIntGen const &x, TryIntGen const &y) {
    return x.val != y.val;
  }
  friend bool operator!=(TryIntGen const &x, int const &y) {
    return x.val != static_cast<Tint>(y);
  }
  //
  friend bool operator>=(TryIntGen const &x, TryIntGen const &y) {
    return x.val >= y.val;
  }
  friend bool operator>=(TryIntGen const &x, int const &y) {
    return x.val >= static_cast<Tint>(y);
  }
  //
  friend bool operator<=(TryIntGen const &x, TryIntGen const &y) {
    return x.val <= y.val;
  }
  friend bool operator<=(TryIntGen const &x, int const &y) {
    return x.val <= static_cast<Tint>(y);
  }
  //
  friend bool operator>(TryIntGen const &x, TryIntGen const &y) {
    return x.val > y.val;
  }
  friend bool operator>(TryIntGen const &x, int const &y) {
    return x.val > static_cast<Tint>(y);
  }
  //
  friend bool operator<(TryIntGen const &x, TryIntGen const &y) {
    return x.val < y.val;
  }
  friend bool operator<(TryIntGen const &x, int const &y) {
    return x.val < static_cast<Tint>(y);
  }
};

// The two concrete flavours, and the backward-compatible default.
using TryCarryInt64 = TryIntGen<try_carry_policy>;
using TrySimdInt64 = TryIntGen<try_simd_policy>;
using TryInt64 = TryCarryInt64;

// Trait: recognizes any deferred-overflow value type (a TryIntGen or a
// Rational over one), so the single terminate_in_arithmetic_error and the
// use_try_int64 exclusion cover both flavours without per-alias boilerplate.
template <typename T> struct is_try_int : std::false_type {};
template <class Pol> struct is_try_int<TryIntGen<Pol>> : std::true_type {};
template <typename T>
struct uses_deferred_overflow : is_try_int<T> {};
template <class Pol>
struct uses_deferred_overflow<Rational<TryIntGen<Pol>>> : std::true_type {};

// The deferred overflow check. Trivial for every type; the deferred-checking
// types check the is_correct flag, reset it, and throw TryIntException so the
// caller can fall back to an exact type. Call it out of the tight loop but not
// completely out of scope (typically once per outer loop iteration) so that a
// computation on wrapped values ends prematurely.
template <typename T> inline void terminate_in_arithmetic_error() {
  if constexpr (uses_deferred_overflow<T>::value) {
    if (!is_correct) {
      is_correct = true;
      throw TryIntException{1};
    }
  }
}

// Conversion from an implementation-of-Z type T (e.g. mpz_class) into the
// target try-type Ttry, throwing TryIntException when the value does not fit
// into int64_t. Some TYPE_CONVERSION into int64_t truncate silently (e.g.
// mpz_class via get_si), so the conversion is verified by a round trip. The
// scratch holds the round-trip value: passing it in from the caller keeps its
// allocation (e.g. the mpz limb) out of entrywise conversion loops. Number
// types with a cheap exact fit test provide a more specialized overload next
// to their own definition (e.g. mpz_class in NumberTheoryGmp.h), which
// overload resolution prefers to this generic version. Ttry defaults to
// TryInt64 so existing call sites need no change.
template <typename Ttry = TryInt64, typename T>
Ttry ConvertToTryInt64(T const &val, T &scratch) {
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
  return Ttry(v);
}

template <typename Ttry = TryInt64, typename T>
Ttry ConvertToTryInt64(T const &val) {
  T scratch;
  return ConvertToTryInt64<Ttry>(val, scratch);
}

template <typename T, class Pol>
T ConvertFromTryInt64(TryIntGen<Pol> const &val) {
  int64_t const &v = val.get_const_val();
  T ret;
  TYPE_CONVERSION(stc<int64_t>{v}, ret);
  return ret;
}

template <class Pol>
inline void TYPE_CONVERSION(stc<int64_t> const &a1, TryIntGen<Pol> &a2) {
  a2 = TryIntGen<Pol>(a1.val);
}

template <class Pol>
inline void TYPE_CONVERSION(stc<TryIntGen<Pol>> const &a1, int64_t &a2) {
  a2 = a1.val.get_const_val();
}

// Traits: same profile as SafeInt64, for both flavours.

template <class Pol> struct is_implementation_of_Z<TryIntGen<Pol>> {
  static const bool value = true;
};

template <class Pol> struct is_implementation_of_Q<TryIntGen<Pol>> {
  static const bool value = false;
};

template <class Pol> struct is_euclidean_domain<TryIntGen<Pol>> {
  static const bool value = true;
};

template <class Pol> struct is_ring_field<TryIntGen<Pol>> {
  static const bool value = false;
};

template <class Pol> struct is_totally_ordered<TryIntGen<Pol>> {
  static const bool value = true;
};

// Exact provided the terminate_in_arithmetic_error discipline is followed: a
// wrapped value is always reported before a result is returned.
template <class Pol> struct is_exact_arithmetic<TryIntGen<Pol>> {
  static const bool value = true;
};

// Native-width integer with no heap allocation: fused form, nothing for a
// scratch to save (see is_fma_prefered).
template <class Pol> struct is_fma_prefered<TryIntGen<Pol>> {
  static const bool value = true;
};

template <class Pol> struct underlying_ring<TryIntGen<Pol>> {
  typedef TryIntGen<Pol> ring_type;
};

// Whether an algorithm should attempt its computation over a try-type first:
// convert the input (throwing TryIntException when an entry does not fit), run
// with the terminate_in_arithmetic_error discipline, and fall back to T on
// TryIntException. The big-integer implementations of Z qualify. The native
// integer types are excluded (already at machine speed, and their callers
// expect the wrapping behavior), and so are the try-types themselves. Types
// whose is_implementation_of_Z is not specialized resolve to false through the
// SFINAE default.
template <typename T, typename = void> struct use_try_int64 {
  static const bool value = false;
};

template <typename T>
struct use_try_int64<T, std::enable_if_t<is_implementation_of_Z<T>::value>> {
  static const bool value =
      !std::is_integral_v<T> && !is_try_int<T>::value;
};

// The Euclidean domain operations, with the same deferred guards as the
// divisions: after an overflow the operands can be wrapped garbage, so a zero
// or overflowing divisor flags is_correct instead of trapping. This is what
// GcdPair / GenericGcd route through, making the gcd-based content reduction
// of the polyhedral kernels available over the try-types.
template <class Pol>
inline void QUO_INT(stc<TryIntGen<Pol>> const &a, stc<TryIntGen<Pol>> const &b,
                    TryIntGen<Pol> &q) {
  using Tint = typename TryIntGen<Pol>::Tint;
  const Tint &a_int = a.val.get_const_val();
  const Tint &b_int = b.val.get_const_val();
  if (b_int == 0 ||
      (b_int == -1 && a_int == std::numeric_limits<int64_t>::min())) {
    is_correct = false;
    q = TryIntGen<Pol>(0);
    return;
  }
  q = TryIntGen<Pol>(QuoInt_C_integer<Tint>(a_int, b_int));
}

template <class Pol>
inline void ResInt_Kernel(TryIntGen<Pol> const &a, TryIntGen<Pol> const &b,
                          TryIntGen<Pol> &res) {
  using Tint = typename TryIntGen<Pol>::Tint;
  const Tint &a_int = a.get_const_val();
  const Tint &b_int = b.get_const_val();
  if (b_int == 0 ||
      (b_int == -1 && a_int == std::numeric_limits<int64_t>::min())) {
    is_correct = false;
    res = TryIntGen<Pol>(0);
    return;
  }
  res.get_val() = ResInt_C_integer<Tint>(a_int, b_int);
}

template <class Pol> inline size_t get_bit(TryIntGen<Pol> const &x) {
  int64_t const &val = x.get_const_val();
  return get_bit(val);
}

// clang-format off
#endif  // SRC_NUMBER_NUMBERTHEORYTRYINT_H_
// clang-format on
