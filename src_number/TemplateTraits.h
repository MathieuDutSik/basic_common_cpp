// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#ifndef SRC_NUMBER_TEMPLATETRAITS_H_
#define SRC_NUMBER_TEMPLATETRAITS_H_

// Trait for specific types.
// the std::is_same<T,mpz_class> is not adequate because it requires the type
// mpz_class to be known in the scope.

#include "PivotCost.h"
#include <cstdint>

template <typename T> struct has_reduction_subset_solver {
  static const bool value = false;
};

template <typename T> struct is_mpq_class {
  static const bool value = false;
};

template <typename T> struct is_mpz_class {
  static const bool value = false;
};

template <typename T> struct is_boost_cpp_int {
  static const bool value = false;
};

template <typename T> struct is_boost_cpp_rational {
  static const bool value = false;
};

template <typename T> struct is_boost_mpz_int {
  static const bool value = false;
};

template <typename T> struct is_boost_mpq_rational {
  static const bool value = false;
};

// is_euclidean_domain

template <typename T> struct is_euclidean_domain {
  static const bool value = false;
};

template <> struct is_euclidean_domain<int64_t> {
  static const bool value = true;
};

template <> struct is_euclidean_domain<int32_t> {
  static const bool value = true;
};

template <> struct is_euclidean_domain<int16_t> {
  static const bool value = true;
};

template <> struct is_euclidean_domain<int8_t> {
  static const bool value = true;
};

// overlying_field block

template <typename T> struct overlying_field {};

template <> struct overlying_field<double> {
  using field_type = double;
};

template <> struct overlying_field<float> {
  using field_type = float;
};

// underlying_ring

template <typename T> struct underlying_ring {};

template <> struct underlying_ring<int64_t> {
  using ring_type = int64_t;
};

template <> struct underlying_ring<int32_t> {
  using ring_type = int32_t;
};

template <> struct underlying_ring<int16_t> {
  using ring_type = int16_t;
};

template <> struct underlying_ring<int8_t> {
  using ring_type = int8_t;
};

// Trait definition for subset of integers

template <typename T> struct is_implementation_of_Z {};

template <> struct is_implementation_of_Z<double> {
  static const bool value = false;
};

template <> struct is_implementation_of_Z<float> {
  static const bool value = false;
};

template <> struct is_implementation_of_Z<int64_t> {
  static const bool value = true;
};

template <> struct is_implementation_of_Z<int32_t> {
  static const bool value = true;
};

template <> struct is_implementation_of_Z<int16_t> {
  static const bool value = true;
};

template <> struct is_implementation_of_Z<int8_t> {
  static const bool value = true;
};

// Whether DeterminantMat should compute the determinant with the Bareiss
// fraction-free algorithm instead of classical Gaussian elimination. Bareiss
// requires an integral domain with EXACT division and controls intermediate
// operand growth, which makes it faster than Gaussian elimination for exact,
// non-trivial arithmetic (integers, rationals, number fields -- benchmarked at
// ~2x for mpq_class and up to ~4x for QuadField). It must be left OFF for:
//   --- floating point, where pivoting for numerical stability matters;
//   --- division-free rings such as jets, which carry zero divisors, so the
//       previous pivot may be non-invertible and the Bareiss division invalid.
// The default follows is_implementation_of_Z (the integer rings); exact fields
// opt in explicitly in their own headers.
template <typename T> struct use_bareiss_for_determinants {
  static const bool value = is_implementation_of_Z<T>::value;
};

// Trait definition for subset of rationals

template <typename T> struct is_implementation_of_Q {};

template <> struct is_implementation_of_Q<double> {
  static const bool value = false;
};

template <> struct is_implementation_of_Q<float> {
  static const bool value = false;
};

template <> struct is_implementation_of_Q<int64_t> {
  static const bool value = false;
};

template <> struct is_implementation_of_Q<int32_t> {
  static const bool value = false;
};

template <> struct is_implementation_of_Q<int16_t> {
  static const bool value = false;
};

template <> struct is_implementation_of_Q<int8_t> {
  static const bool value = false;
};

// Trait definition for exactness

template <typename T> struct is_exact_arithmetic {};

template <> struct is_exact_arithmetic<double> {
  static const bool value = false;
};

template <> struct is_exact_arithmetic<float> {
  static const bool value = false;
};

// Trait for the preferred multiply-accumulate ("acc += a * b") form.
//
// value == true  (the default): the direct fused form  acc += a * b  is at
//   least as fast as materializing the product into a temporary first. True for
//   native types, for boost::multiprecision::mpz_int (whose expression
//   templates fuse += product), and for the number types whose
//   operator+=(product-proxy) accumulates in place (the jet / RealField
//   expression templates). The generic matrix kernels then keep writing
//   acc += a * b.
//
// value == false: a reused scratch is measurably faster,
//   prod = a * b; acc += prod;
//   because acc += a * b would allocate and free a fresh temporary on every
//   evaluation (gmpxx / boost.cpp: mpz_class, mpq_class, cpp_int, cpp_rational,
//   mpq_rational), or because operator=(product) is cheaper than
//   operator+=(product) for the type (QuadField). Both branches compute the
//   same value; the trait only selects the faster implementation.
template <typename T> struct is_fma_prefered {
  static const bool value = true;
};

// Empty placeholder for a reuse-scratch that a code path does not need. When
// is_fma_prefered<T> is true a kernel declares its scratch as
//   std::conditional_t<is_fma_prefered<T>::value, empty_scratch, T> scratch;
// so that for the fused-preferring types NO unused T (e.g. an mpq_class or a
// jet) is constructed at all -- the scratch collapses to this empty object,
// which the compiler need not even lay out.
struct empty_scratch {};

// Native types: the product stays in a register, so the direct/fused form is
// best (a scratch has nothing to save). These match the default but are stated
// explicitly so every numerical type carries a deliberate choice.
template <> struct is_fma_prefered<int8_t> {
  static const bool value = true;
};
template <> struct is_fma_prefered<int16_t> {
  static const bool value = true;
};
template <> struct is_fma_prefered<int32_t> {
  static const bool value = true;
};
template <> struct is_fma_prefered<int64_t> {
  static const bool value = true;
};
template <> struct is_fma_prefered<double> {
  static const bool value = true;
};
template <> struct is_fma_prefered<float> {
  static const bool value = true;
};

// Trait definition for fields

template <typename T> struct is_ring_field {};

template <> struct is_ring_field<int64_t> {
  static const bool value = false;
};

template <> struct is_ring_field<int32_t> {
  static const bool value = false;
};

template <> struct is_ring_field<int16_t> {
  static const bool value = false;
};

template <> struct is_ring_field<int8_t> {
  static const bool value = false;
};

template <> struct is_ring_field<double> {
  static const bool value = true;
};

template <> struct is_ring_field<float> {
  static const bool value = true;
};

// Trait selecting a division-free determinant algorithm. For fields and integral
// domains the elimination-based determinant (DeterminantMatKernel) is exact and
// fast. For rings with zero divisors -- in particular the truncated jet ring
// T[t]/(t^{N+1}), where a rationally degenerate quantity has a zero constant
// term and is therefore not invertible -- elimination can be forced to divide by
// a zero divisor. Types that set this to true are routed through a division-free
// determinant instead. Default: false.
template <typename T> struct determinant_division_free {
  static const bool value = false;
};

// The value of a scalar at its degeneracy / base point. A generic scalar is its
// own constant term (identity); the jet specialization (in jet_number.h) returns
// the t = 0 coefficient. This is the bridge that lets generic code (e.g. the
// division-free determinant dispatch) build the constant-term matrix.
template <typename T> T const &constant_term(T const &x) { return x; }

// Trait of totally ordered set

template <typename T> struct is_totally_ordered {
  static const bool value = false;
};

template <> struct is_totally_ordered<int64_t> {
  static const bool value = true;
};

template <> struct is_totally_ordered<int32_t> {
  static const bool value = true;
};

template <> struct is_totally_ordered<int16_t> {
  static const bool value = true;
};

template <> struct is_totally_ordered<int8_t> {
  static const bool value = true;
};

template <> struct is_totally_ordered<double> {
  static const bool value = true;
};

template <> struct is_totally_ordered<float> {
  static const bool value = true;
};

// is floating arithmetic

template <typename T> struct is_float_arithmetic {
  static const bool value = false;
};

template <> struct is_float_arithmetic<float> {
  static const bool value = true;
};

template <> struct is_float_arithmetic<double> {
  static const bool value = true;
};

// Trait definition for underlying ring

template <typename T> struct underlying_totally_ordered_ring {};

template <> struct underlying_totally_ordered_ring<int64_t> {
  using real_type = int64_t;
};

template <> struct underlying_totally_ordered_ring<int32_t> {
  using real_type = int32_t;
};

template <> struct underlying_totally_ordered_ring<int16_t> {
  using real_type = int16_t;
};

template <> struct underlying_totally_ordered_ring<int8_t> {
  using real_type = int8_t;
};

// clang-format off
#endif  // SRC_NUMBER_TEMPLATETRAITS_H_
// clang-format on
