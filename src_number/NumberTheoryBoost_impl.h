// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
/*
  The wrappers of a pair of boost::multiprecision types, an integer type and
  its rational type, shared by NumberTheoryBoostCppInt.h (cpp_int,
  cpp_rational) and NumberTheoryBoostGmpInt.h (mpz_int, mpq_rational). Those
  two headers were the same code but for the type names; each now defines
  the names below and includes this file, then adds what is specific to its
  backend.

  No include guard: it is included once per pair of types. The macros are
  undefined at the end.
  --- NT_BOOST_INT, NT_BOOST_RAT: the two types.
  --- NT_BOOST_IS_INT, NT_BOOST_IS_RAT: the traits naming them.
  --- NT_BOOST_INT_FMA_PREFERED: is_fma_prefered of the integer type.
  --- NT_BOOST_INT_TO_SMALL, NT_BOOST_CEIL_RAT, NT_BOOST_FLOOR_RAT: the
      names of three helpers, kept per type.
 */
#if !defined(NT_BOOST_INT) || !defined(NT_BOOST_RAT) ||                       \
    !defined(NT_BOOST_IS_INT) || !defined(NT_BOOST_IS_RAT) ||                 \
    !defined(NT_BOOST_INT_FMA_PREFERED) || !defined(NT_BOOST_INT_TO_SMALL) ||  \
    !defined(NT_BOOST_CEIL_RAT) || !defined(NT_BOOST_FLOOR_RAT)
#error "NumberTheoryBoost_impl.h needs its NT_BOOST_* macros"
#endif

template <> struct NT_BOOST_IS_INT<NT_BOOST_INT> {
  static const bool value = true;
};

template <> struct NT_BOOST_IS_RAT<NT_BOOST_RAT> {
  static const bool value = true;
};

// hash

namespace std {
template <> struct hash<NT_BOOST_INT> {
  std::size_t operator()(const NT_BOOST_INT &val) const {
    std::stringstream s;
    s << val;
    std::string converted(s.str());
    return std::hash<std::string>()(converted);
  }
};
template <> struct hash<NT_BOOST_RAT> {
  std::size_t operator()(const NT_BOOST_RAT &val) const {
    std::stringstream s;
    s << val;
    std::string converted(s.str());
    return std::hash<std::string>()(converted);
  }
};
// clang-format off
}  // namespace std
// clang-format on


template <> struct is_euclidean_domain<NT_BOOST_INT> {
  static const bool value = true;
};
template <> struct is_euclidean_domain<NT_BOOST_RAT> {
  static const bool value = true;
};

template <> struct is_exact_arithmetic<NT_BOOST_INT> {
  static const bool value = true;
};
template <> struct is_exact_arithmetic<NT_BOOST_RAT> {
  static const bool value = true;
};

template <> struct is_implementation_of_Z<NT_BOOST_INT> {
  static const bool value = true;
};
template <> struct is_implementation_of_Z<NT_BOOST_RAT> {
  static const bool value = false;
};

// FMA form (see is_fma_prefered): per backend, see NT_BOOST_INT_FMA_PREFERED.
// The rational type materializes a temporary for `a*b`, so the
// reused-scratch form wins.
template <> struct is_fma_prefered<NT_BOOST_INT> {
  static const bool value = NT_BOOST_INT_FMA_PREFERED;
};
template <> struct is_fma_prefered<NT_BOOST_RAT> {
  static const bool value = false;
};

template <> struct is_implementation_of_Q<NT_BOOST_INT> {
  static const bool value = false;
};
template <> struct is_implementation_of_Q<NT_BOOST_RAT> {
  static const bool value = true;
};

template <> struct is_ring_field<NT_BOOST_INT> {
  static const bool value = false;
};
template <> struct is_ring_field<NT_BOOST_RAT> {
  static const bool value = true;
};
// Exact fraction field: opt into Bareiss (see use_bareiss_for_determinants).
template <>
struct use_bareiss_for_determinants<NT_BOOST_RAT> {
  static const bool value = true;
};

template <> struct is_totally_ordered<NT_BOOST_INT> {
  static const bool value = true;
};

template <> struct is_totally_ordered<NT_BOOST_RAT> {
  static const bool value = true;
};

template <> struct underlying_ring<NT_BOOST_INT> {
  typedef NT_BOOST_INT ring_type;
};
template <> struct underlying_ring<NT_BOOST_RAT> {
  typedef NT_BOOST_INT ring_type;
};

template <> struct underlying_z_ring<NT_BOOST_INT> {
  typedef NT_BOOST_INT ring_type;
};

template <> struct underlying_z_ring<NT_BOOST_RAT> {
  typedef NT_BOOST_INT ring_type;
};

// the integer type is a ring and contains no Q, so it gets no underlying_q_field.
template <> struct underlying_q_field<NT_BOOST_RAT> {
  typedef NT_BOOST_RAT field_type;
};

template <> struct overlying_field<NT_BOOST_INT> {
  typedef NT_BOOST_RAT field_type;
};
template <> struct overlying_field<NT_BOOST_RAT> {
  typedef NT_BOOST_RAT field_type;
};

template <>
struct underlying_totally_ordered_ring<NT_BOOST_INT> {
  typedef NT_BOOST_INT real_type;
};
template <>
struct underlying_totally_ordered_ring<NT_BOOST_RAT> {
  typedef NT_BOOST_RAT real_type;
};

inline NT_BOOST_INT
CanonicalizationUnit(NT_BOOST_INT const &eVal) {
  if (eVal < 0)
    return -1;
  return 1;
}
inline NT_BOOST_RAT
CanonicalizationUnit(NT_BOOST_RAT const &eVal) {
  if (eVal < 0)
    return -1;
  return 1;
}

namespace boost::multiprecision {
inline void ResInt_Kernel(NT_BOOST_INT const &a,
                          NT_BOOST_INT const &b,
                          NT_BOOST_INT &res) {
  using T = NT_BOOST_INT;
  T q = a / b;
  if (a < 0 && b * q != a) {
    if (b > 0)
      q--;
    else
      q++;
  }
  res = a - q * b;
}
}  // namespace boost::multiprecision

void QUO_INT(stc<NT_BOOST_INT> const &a,
             stc<NT_BOOST_INT> const &b,
             NT_BOOST_INT &q) {
  q = a.val / b.val;
  if (a.val < 0 && b.val * q != a.val) {
    if (b.val > 0)
      q--;
    else
      q++;
  }
}

inline std::pair<NT_BOOST_RAT,
                 NT_BOOST_RAT>
ResQuoInt_kernel(NT_BOOST_RAT const &a,
                 NT_BOOST_RAT const &b) {
  // a = a_n / a_d
  // b = b_n / b_d
  // a = res + q * b  with 0 <= res < |b|
  // equivalent to
  // a_n / a_d = res + q * (b_n / b_d)
  // equivalent to
  // a_n * b_d = res * a_d * b_d + (q * a_d) * b_n
  using Tf = NT_BOOST_RAT;
  using T = NT_BOOST_INT;
  T a_n = numerator(a);
  T b_n = numerator(b);
  T a_d = denominator(a);
  T b_d = denominator(b);
  T a1 = a_n * b_d;
  T b1 = a_d * b_n;
  T q = a1 / b1;
  Tf q_f = q;
  Tf res = a - q_f * b;
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
          std::cerr << "Error in ResQuoInt_kernel for "
                       "NT_BOOST_RAT\n";
          throw TerminalException{1};
        }
        return {res, q};
      }
    }
  }
}
namespace boost::multiprecision {
inline void ResInt_Kernel(NT_BOOST_RAT const &a,
                          NT_BOOST_RAT const &b,
                          NT_BOOST_RAT &res) {
  res = ResQuoInt_kernel(a, b).first;
}
}  // namespace boost::multiprecision

void QUO_INT(stc<NT_BOOST_RAT> const &a,
             stc<NT_BOOST_RAT> const &b,
             NT_BOOST_RAT &q) {
  q = ResQuoInt_kernel(a.val, b.val).second;
}

#include "QuoIntFcts.h"

inline bool IsInteger(NT_BOOST_RAT const &x) {
  NT_BOOST_INT one = 1;
  NT_BOOST_INT eDen = denominator(x);
  return eDen == one;
}

//

inline NT_BOOST_INT
GetNumerator(NT_BOOST_INT const &x) {
  return x;
}

inline NT_BOOST_RAT
GetNumerator(NT_BOOST_RAT const &x) {
  NT_BOOST_INT eNum = numerator(x);
  NT_BOOST_RAT eNum_q = eNum;
  return eNum_q;
}

//

inline NT_BOOST_INT
GetDenominator([[maybe_unused]] NT_BOOST_INT const &x) {
  return 1;
}

inline NT_BOOST_RAT
GetDenominator(NT_BOOST_RAT const &x) {
  NT_BOOST_INT eDen = denominator(x);
  NT_BOOST_RAT eDen_q = eDen;
  return eDen_q;
}

//

inline NT_BOOST_INT
GetNumerator_z(NT_BOOST_RAT const &x) {
  return numerator(x);
}

inline NT_BOOST_INT
GetNumerator_z(NT_BOOST_INT const &x) {
  return x;
}

//

inline NT_BOOST_INT
GetDenominator_z(NT_BOOST_RAT const &x) {
  return denominator(x);
}

inline NT_BOOST_INT
GetDenominator_z([[maybe_unused]] NT_BOOST_INT const &x) {
  return 1;
}

//

inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            double &a2) {
  a2 = a1.val.template convert_to<double>();
}

// double as input.
// This mirrors the conversion in NumberTheoryGmp.h: the double value is
// truncated to int64_t before being assigned to the big integer. This
// only handles doubles whose truncated value fits in int64_t; conversion
// from doubles outside that range has inherent issues that must be
// addressed by the caller.
inline void TYPE_CONVERSION(stc<double> const &a1,
                            NT_BOOST_INT &a2) {
  a2 = static_cast<int64_t>(a1.val);
}
template <typename T>
  requires (std::is_same_v<T, int>
            && !std::is_same_v<int, int8_t>
            && !std::is_same_v<int, int16_t>
            && !std::is_same_v<int, int32_t>
            && !std::is_same_v<int, int64_t>)
inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            T &a2) {
  NT_BOOST_INT_TO_SMALL(a1.val, a2);
}
template <typename T>
  requires (std::is_same_v<T, long>
            && !std::is_same_v<long, int8_t>
            && !std::is_same_v<long, int16_t>
            && !std::is_same_v<long, int32_t>
            && !std::is_same_v<long, int64_t>)
inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            T &a2) {
  NT_BOOST_INT_TO_SMALL(a1.val, a2);
}
template <typename T>
  requires (std::is_same_v<T, int>
            && !std::is_same_v<int, int8_t>
            && !std::is_same_v<int, int16_t>
            && !std::is_same_v<int, int32_t>
            && !std::is_same_v<int, int64_t>)
inline void TYPE_CONVERSION(stc<T> const &a1,
                            NT_BOOST_INT &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            NT_BOOST_INT &a2) {
  a2 = a1.val;
}

inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            double &a2) {
  a2 = a1.val.template convert_to<double>();
}

inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            NT_BOOST_INT &a2) {
  if (!IsInteger(a1.val)) {
    std::string str = std::format("a1={} is not an integer", a1.val);
    throw ConversionException{str};
  }
  a2 = numerator(a1.val);
}
template <typename T>
  requires (std::is_same_v<T, int>
            && !std::is_same_v<int, int8_t>
            && !std::is_same_v<int, int16_t>
            && !std::is_same_v<int, int32_t>
            && !std::is_same_v<int, int64_t>)
inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            T &a2) {
  NT_BOOST_INT a1_z;
  TYPE_CONVERSION(a1, a1_z);
  stc<NT_BOOST_INT> stc_a1_z{a1_z};
  TYPE_CONVERSION(stc_a1_z, a2);
}
template <typename T>
  requires (std::is_same_v<T, long>
            && !std::is_same_v<long, int8_t>
            && !std::is_same_v<long, int16_t>
            && !std::is_same_v<long, int32_t>
            && !std::is_same_v<long, int64_t>)
inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            T &a2) {
  NT_BOOST_INT a1_z;
  TYPE_CONVERSION(a1, a1_z);
  stc<NT_BOOST_INT> stc_a1_z{a1_z};
  TYPE_CONVERSION(stc_a1_z, a2);
}
template <typename T>
  requires (std::is_same_v<T, int>
            && !std::is_same_v<int, int8_t>
            && !std::is_same_v<int, int16_t>
            && !std::is_same_v<int, int32_t>
            && !std::is_same_v<int, int64_t>)
inline void TYPE_CONVERSION(stc<T> const &a1,
                            NT_BOOST_RAT &a2) {
  a2 = a1.val;
}
template <typename T>
  requires (std::is_same_v<T, long>
            && !std::is_same_v<long, int8_t>
            && !std::is_same_v<long, int16_t>
            && !std::is_same_v<long, int32_t>
            && !std::is_same_v<long, int64_t>)
inline void TYPE_CONVERSION(stc<T> const &a1,
                            NT_BOOST_RAT &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            NT_BOOST_RAT &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            NT_BOOST_RAT &a2) {
  a2 = a1.val;
}

// int8_t as input

inline void TYPE_CONVERSION(stc<int8_t> const &a1,
                            NT_BOOST_INT &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<int8_t> const &a1,
                            NT_BOOST_RAT &a2) {
  a2 = a1.val;
}

// uint8_t as input

inline void TYPE_CONVERSION(stc<uint8_t> const &a1,
                            NT_BOOST_INT &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<uint8_t> const &a1,
                            NT_BOOST_RAT &a2) {
  a2 = a1.val;
}

// int16_t as input

inline void TYPE_CONVERSION(stc<int16_t> const &a1,
                            NT_BOOST_INT &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<int16_t> const &a1,
                            NT_BOOST_RAT &a2) {
  a2 = a1.val;
}

// uint16_t as input

inline void TYPE_CONVERSION(stc<uint16_t> const &a1,
                            NT_BOOST_INT &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<uint16_t> const &a1,
                            NT_BOOST_RAT &a2) {
  a2 = a1.val;
}

// uint32_t as input

inline void TYPE_CONVERSION(stc<uint32_t> const &a1,
                            NT_BOOST_INT &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<uint32_t> const &a1,
                            NT_BOOST_RAT &a2) {
  a2 = a1.val;
}

// int32_t as input.
inline void TYPE_CONVERSION(stc<int32_t> const &a1,
                            NT_BOOST_INT &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<int32_t> const &a1,
                            NT_BOOST_RAT &a2) {
  a2 = a1.val;
}

// int64_t as input.
inline void TYPE_CONVERSION(stc<int64_t> const &a1,
                            NT_BOOST_INT &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<int64_t> const &a1,
                            NT_BOOST_RAT &a2) {
  a2 = a1.val;
}

// long as input: only enabled when long is distinct from every fixed-width
// integer type for which an overload is already defined above.
template <typename T>
  requires (std::is_same_v<T, long>
            && !std::is_same_v<long, int8_t>
            && !std::is_same_v<long, int16_t>
            && !std::is_same_v<long, int32_t>
            && !std::is_same_v<long, int64_t>)
inline void TYPE_CONVERSION(stc<T> const &a1,
                            NT_BOOST_INT &a2) {
  a2 = a1.val;
}

// uint64_t as input.
// No `unsigned long` overloads are defined elsewhere, so no platform guard
// is needed: on Linux uint64_t == unsigned long, and on Apple they are
// distinct types but unsigned long is not used as an overload elsewhere.
inline void TYPE_CONVERSION(stc<uint64_t> const &a1,
                            NT_BOOST_INT &a2) {
  a2 = a1.val;
}
inline void TYPE_CONVERSION(stc<uint64_t> const &a1,
                            NT_BOOST_RAT &a2) {
  a2 = a1.val;
}

// size_t as input.
// On platforms where size_t coincides with uint64_t (e.g. Linux x86_64)
// or with uint32_t (typical 32-bit platforms), the existing non-template
// overloads above already cover the call and adding plain overloads here
// would cause redefinition errors. We therefore guard these with a
// `requires` clause that excludes the size_t == uint64_t and size_t ==
// uint32_t cases; on those platforms the templates produce no candidates
// and the non-template overloads win. On platforms where size_t is a
// distinct type (e.g. Apple, where size_t is `unsigned long` while
// uint64_t is `unsigned long long`), the templates supply the needed
// overload.
template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1,
                            NT_BOOST_INT &a2) {
  a2 = a1.val;
}
template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<T> const &a1,
                            NT_BOOST_RAT &a2) {
  a2 = a1.val;
}

// the integer type as output target for the small integer types.
// The boost convert_to<>() does not check overflow, so we perform an
// explicit range check against the destination type before converting
// and throw ConversionException if the value does not fit.

template <typename Tout>
inline void NT_BOOST_INT_TO_SMALL(NT_BOOST_INT const &val,
                                     Tout &out) {
  if (val < std::numeric_limits<Tout>::min() ||
      val > std::numeric_limits<Tout>::max()) {
    std::string str = "value=" + val.str() +
                      " does not fit in the destination integer type";
    throw ConversionException{str};
  }
  out = val.template convert_to<Tout>();
}

inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            int8_t &a2) {
  NT_BOOST_INT_TO_SMALL(a1.val, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            uint8_t &a2) {
  NT_BOOST_INT_TO_SMALL(a1.val, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            int16_t &a2) {
  NT_BOOST_INT_TO_SMALL(a1.val, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            uint16_t &a2) {
  NT_BOOST_INT_TO_SMALL(a1.val, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            int32_t &a2) {
  NT_BOOST_INT_TO_SMALL(a1.val, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            uint32_t &a2) {
  NT_BOOST_INT_TO_SMALL(a1.val, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            int64_t &a2) {
  NT_BOOST_INT_TO_SMALL(a1.val, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            uint64_t &a2) {
  NT_BOOST_INT_TO_SMALL(a1.val, a2);
}
// size_t output: only enabled when size_t is a distinct type from
// uint64_t and uint32_t (see the size_t input section for rationale).
template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<NT_BOOST_INT> const &a1,
                            T &a2) {
  NT_BOOST_INT_TO_SMALL(a1.val, a2);
}

// the rational type as output target for the small integer types

inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            int8_t &a2) {
  NT_BOOST_INT a1_z;
  TYPE_CONVERSION(a1, a1_z);
  stc<NT_BOOST_INT> stc_a1_z{a1_z};
  TYPE_CONVERSION(stc_a1_z, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            uint8_t &a2) {
  NT_BOOST_INT a1_z;
  TYPE_CONVERSION(a1, a1_z);
  stc<NT_BOOST_INT> stc_a1_z{a1_z};
  TYPE_CONVERSION(stc_a1_z, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            int16_t &a2) {
  NT_BOOST_INT a1_z;
  TYPE_CONVERSION(a1, a1_z);
  stc<NT_BOOST_INT> stc_a1_z{a1_z};
  TYPE_CONVERSION(stc_a1_z, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            uint16_t &a2) {
  NT_BOOST_INT a1_z;
  TYPE_CONVERSION(a1, a1_z);
  stc<NT_BOOST_INT> stc_a1_z{a1_z};
  TYPE_CONVERSION(stc_a1_z, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            int32_t &a2) {
  NT_BOOST_INT a1_z;
  TYPE_CONVERSION(a1, a1_z);
  stc<NT_BOOST_INT> stc_a1_z{a1_z};
  TYPE_CONVERSION(stc_a1_z, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            uint32_t &a2) {
  NT_BOOST_INT a1_z;
  TYPE_CONVERSION(a1, a1_z);
  stc<NT_BOOST_INT> stc_a1_z{a1_z};
  TYPE_CONVERSION(stc_a1_z, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            int64_t &a2) {
  NT_BOOST_INT a1_z;
  TYPE_CONVERSION(a1, a1_z);
  stc<NT_BOOST_INT> stc_a1_z{a1_z};
  TYPE_CONVERSION(stc_a1_z, a2);
}
inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            uint64_t &a2) {
  NT_BOOST_INT a1_z;
  TYPE_CONVERSION(a1, a1_z);
  stc<NT_BOOST_INT> stc_a1_z{a1_z};
  TYPE_CONVERSION(stc_a1_z, a2);
}
// size_t output: only enabled when size_t is a distinct type from
// uint64_t and uint32_t (see the size_t input section for rationale).
template <typename T>
  requires (std::is_same_v<T, size_t>
            && !std::is_same_v<size_t, uint64_t>
            && !std::is_same_v<size_t, uint32_t>)
inline void TYPE_CONVERSION(stc<NT_BOOST_RAT> const &a1,
                            T &a2) {
  NT_BOOST_INT a1_z;
  TYPE_CONVERSION(a1, a1_z);
  stc<NT_BOOST_INT> stc_a1_z{a1_z};
  TYPE_CONVERSION(stc_a1_z, a2);
}

inline void
ScalingInteger_Kernel(stc<NT_BOOST_RAT> const &x,
                      NT_BOOST_INT &x_ret) {
  x_ret = denominator(x.val);
}

inline void ScalingInteger_Kernel(
    [[maybe_unused]] stc<NT_BOOST_INT> const &x,
    NT_BOOST_INT &x_ret) {
  x_ret = 1;
}

inline NT_BOOST_RAT
FractionalPart(NT_BOOST_RAT const &x) {
  using T = NT_BOOST_INT;
  using Tf = NT_BOOST_RAT;
  T x_n = numerator(x);
  T x_d = denominator(x);
  T res;
  ResInt_Kernel(x_n, x_d, res);
  Tf res_f = res;
  Tf x_df = x_d;
  Tf ret = res_f / x_df;
  return ret;
}

inline NT_BOOST_RAT
NT_BOOST_FLOOR_RAT(NT_BOOST_RAT const &x) {
  NT_BOOST_RAT eFrac = FractionalPart(x);
  return x - eFrac;
}

inline NT_BOOST_RAT
NT_BOOST_CEIL_RAT(NT_BOOST_RAT const &x) {
  NT_BOOST_RAT eFrac = FractionalPart(x);
  if (eFrac == 0)
    return x;
  return 1 + x - eFrac;
}

inline void FloorInteger(NT_BOOST_RAT const &xI,
                         NT_BOOST_RAT &xO) {
  xO = NT_BOOST_FLOOR_RAT(xI);
}
inline void FloorInteger(NT_BOOST_RAT const &xI,
                         NT_BOOST_INT &xO) {
  xO = numerator(NT_BOOST_FLOOR_RAT(xI));
}
inline void FloorInteger(NT_BOOST_RAT const &xI,
                         int &xO) {
  NT_BOOST_INT val = numerator(NT_BOOST_FLOOR_RAT(xI));
  xO = val.template convert_to<int>();
}
inline void FloorInteger(NT_BOOST_RAT const &xI,
                         long &xO) {
  NT_BOOST_INT val = numerator(NT_BOOST_FLOOR_RAT(xI));
  xO = val.template convert_to<long>();
}

inline void CeilInteger(NT_BOOST_RAT const &xI,
                        NT_BOOST_RAT &xO) {
  xO = NT_BOOST_CEIL_RAT(xI);
}
inline void CeilInteger(NT_BOOST_RAT const &xI,
                        NT_BOOST_INT &xO) {
  xO = numerator(NT_BOOST_CEIL_RAT(xI));
}
inline void CeilInteger(NT_BOOST_RAT const &xI,
                        int &xO) {
  NT_BOOST_INT val = numerator(NT_BOOST_CEIL_RAT(xI));
  xO = val.template convert_to<int>();
}
inline void CeilInteger(NT_BOOST_RAT const &xI,
                        long &xO) {
  NT_BOOST_INT val = numerator(NT_BOOST_CEIL_RAT(xI));
  xO = val.template convert_to<long>();
}

inline NT_BOOST_RAT
NearestInteger_rni(NT_BOOST_RAT const &x) {
  NT_BOOST_RAT eFrac = FractionalPart(x);
  NT_BOOST_RAT eDiff1 = eFrac;
  NT_BOOST_RAT eDiff2 = 1 - eFrac;
  NT_BOOST_RAT RetVal = x - eFrac;
  if (eDiff1 <= eDiff2) {
    return RetVal;
  } else {
    return RetVal + 1;
  }
}
inline void NearestInteger(NT_BOOST_RAT const &xI,
                           NT_BOOST_RAT &xO) {
  xO = NearestInteger_rni(xI);
}
inline void NearestInteger(NT_BOOST_RAT const &xI,
                           NT_BOOST_INT &xO) {
  NT_BOOST_RAT xO_q = NearestInteger_rni(xI);
  xO = numerator(xO_q);
}

inline void set_to_infinity(NT_BOOST_RAT &x) {
  x = std::numeric_limits<uint64_t>::max();
}

inline void set_to_infinity(NT_BOOST_INT &x) {
  x = std::numeric_limits<uint64_t>::max();
}

inline NT_BOOST_RAT
T_NormGen(NT_BOOST_RAT const &x) {
  if (x < 0)
    return -x;
  return x;
}

inline NT_BOOST_INT
T_NormGen(NT_BOOST_INT const &x) {
  if (x < 0)
    return -x;
  return x;
}

bool universal_square_root(NT_BOOST_INT &ret,
                           NT_BOOST_INT const &val) {
  using T = NT_BOOST_INT;
  ret = sqrt(val);
  T eProd = ret * ret;
  return eProd == val;
}

bool universal_square_root(NT_BOOST_RAT &ret,
                           NT_BOOST_RAT const &val) {
  using T = NT_BOOST_INT;
  using Tf = NT_BOOST_RAT;
  T val_n = numerator(val);
  T val_d = denominator(val);
  T ret_n, ret_d;
  if (!universal_square_root(ret_n, val_n))
    return false;
  if (!universal_square_root(ret_d, val_d))
    return false;
  ret = Tf(ret_n) / Tf(ret_d);
  return true;
}

#undef NT_BOOST_INT
#undef NT_BOOST_RAT
#undef NT_BOOST_IS_INT
#undef NT_BOOST_IS_RAT
#undef NT_BOOST_INT_FMA_PREFERED
#undef NT_BOOST_INT_TO_SMALL
#undef NT_BOOST_CEIL_RAT
#undef NT_BOOST_FLOOR_RAT
