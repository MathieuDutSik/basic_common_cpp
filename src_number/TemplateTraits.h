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

// The flint types (NumberTheoryFlint.h, with ENABLE_FLINT_SUPPORT). The
// primary templates live here so that the matrix code can test them
// without the flint headers in scope, exactly like is_mpz_class.

template <typename T> struct is_fmpz_class {
  static const bool value = false;
};

template <typename T> struct is_fmpq_class {
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

// underlying_z_ring and underlying_q_field
//
// The rational integers Z and the rational numbers Q sitting inside T. These
// answer a different question from underlying_ring, and the two only look
// alike because they agree on every type whose scalars are already rational.
//
//   underlying_ring<T>    a ring inside T over which the computation can run
//                         without denominators, for Bareiss and the other
//                         fraction-free paths. It stays as close to T as it
//                         can: for QuadField<mpq_class, d> it is
//                         Z[sqrt(d)] = QuadField<mpz_class, d>, and for a real
//                         algebraic field it is the order Z[x].
//   underlying_z_ring<T>  the rational integers inside T, so mpz_class for
//                         mpq_class and equally for QuadField<mpq_class, d> or
//                         RealField: the algebraic extension is left behind.
//   underlying_q_field<T> the rational numbers inside T, mpq_class in the same
//                         three cases.
//
// The distinction matters wherever the object being computed is a lattice over
// Z rather than a module over the ring of T -- the basis transformation of an
// LLL reduction, the index of a sublattice being factored, the integral lift
// of a permutation of short vectors. All of those are rational whatever field
// the quadratic form takes its values in, and asking underlying_ring for them
// lands in Z[sqrt(d)], which is not even a euclidean domain for most d.
//
// underlying_q_field is left undefined for a ring, which does not contain Q,
// and both are left undefined for a type with no rational scalars at all --
// double and float, the finite field Fp, ThresholdField and jet. That is a
// hard error rather than a silent answer, as for the other traits of this
// file. Code that has to branch can detect the absence:
//
//   template <typename T>
//   concept HasUnderlyingQField = requires {
//     typename underlying_q_field<T>::field_type;
//   };

template <typename T> struct underlying_z_ring {};

template <> struct underlying_z_ring<int64_t> {
  using ring_type = int64_t;
};

template <> struct underlying_z_ring<int32_t> {
  using ring_type = int32_t;
};

template <> struct underlying_z_ring<int16_t> {
  using ring_type = int16_t;
};

template <> struct underlying_z_ring<int8_t> {
  using ring_type = int8_t;
};

template <typename T> struct underlying_q_field {};

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

// Whether Inverse should be computed with the fraction-free LU factorization
// (Zhou & Jeffrey, InverseFractionFreeLU) instead of classical Gauss-Jordan.
// The forward-elimination-then-back-substitution structure does fewer ring
// operations than the Bareiss-Montante Gauss-Jordan and, benchmarked, beats
// classical elimination for the integer rings and for the exact number fields
// (QuadField, RealField), where it wins by ~1.3x (integers) to several times
// (number fields). It is left OFF for:
//   --- mpq_class and the rational fields, where classical is competitive and
//       the fraction-free gain only appears in dimensions we rarely reach;
//   --- floating point, where pivoting for numerical stability matters.
// The default follows is_implementation_of_Z (the integer rings); the exact
// number fields opt in explicitly in their own headers.
template <typename T> struct use_fraction_free_lu {
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
//   same value; the trait only selects the faster implementation. The
//   choices were re-measured with the GMP pool installed and are
//   unchanged by it.
//
// The kernels consume the trait through AddMul / SubMul below rather
// than branching themselves.
template <typename T> struct is_fma_prefered {
  static const bool value = true;
};

// The multiply-accumulate acc += a * b and acc -= a * b in the fastest
// form for the type: the compound expression when it is at least as
// good (is_fma_prefered), a reused thread local scratch otherwise, and
// the native fused calls where they exist (mpz_class, specialized in
// NumberTheoryGmp.h: the gmpxx expression templates do NOT fuse the
// compound form, mpz_addmul / mpz_submul must be called directly).
template <typename T>
inline void AddMul(T &acc, T const &a, T const &b) {
  if constexpr (is_fma_prefered<T>::value) {
    acc += a * b;
  } else {
    static thread_local T scratch;
    scratch = a * b;
    acc += scratch;
  }
}

template <typename T>
inline void SubMul(T &acc, T const &a, T const &b) {
  if constexpr (is_fma_prefered<T>::value) {
    acc -= a * b;
  } else {
    static thread_local T scratch;
    scratch = a * b;
    acc -= scratch;
  }
}

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

// Whether the ring canonicalizes a vector up to a positive scalar on its own,
// through ScalarCanonicalizationVectorRing. A ring that is neither Z-like nor
// a field has no gcd to reduce a content with and no division to normalize
// with, so it otherwise has to go through its overlying field, which costs a
// conversion and a field division per vector.
template <typename T> struct has_ring_canonicalization {
  static const bool value = false;
};

// Trait definition for real algebraic fields.
//
// The primary template lives here, with the other traits, rather than next to
// RealField: NumberTheoryThreshold.h specializes it without depending on
// NumberTheoryRealField.h, and a specialization of an undeclared template is
// an error.
//
// Unlike is_ring_field and the other traits of this file, this one answers
// "is this one of the real algebraic types" and so has to be false for every
// type that never specializes it. An empty primary would not do: the trait is
// read from the requires clauses guarding the conversions out of RealField and
// RealRing, and a missing value member there is a substitution failure, which
// silently drops the conversion instead of selecting it.

template <typename T> struct is_real_algebraic_field {
  static const bool value = false;
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

// Whether ComputeRowHermiteNormalForm_second uses the modulo-D algorithm
// (Domich-Kannan-Trotter, see HermiteNormalFormModD_or_none): all the
// intermediate entries stay bounded by (a multiple of) the lattice
// determinant D, which removes the coefficient explosion of the naive
// elimination. The requirements are an exact euclidean domain with
// ComputePairGcdDot / ResInt / QuoInt, so the types opt in explicitly in
// their own headers (mpz_class, fmpz_class); the bounded types stay out
// because the Bareiss determinant D overflows before the naive elimination
// would.
template <typename T> struct use_hnf_mod_D {
  static const bool value = false;
};

// Whether SmithNormalFormInvariant uses the Kannan-Bachem alternation of
// modulo-D Hermite reductions (SmithNormalFormInvariantModD_or_none).
// Benchmarked OFF for mpz_class: unlike the row-only Hermite elimination,
// the generic Smith kernel reduces from both sides and shows no
// coefficient explosion, and it beats the alternation by ~1.2-2x at every
// size probed (16..96, spreads 10..10^6). The machinery stays available
// for a ring whose generic Smith kernel does blow up: opting in is this
// one trait.
template <typename T> struct use_snf_mod_D {
  static const bool value = false;
};

// Whether SmithNormalFormInvariant runs the sparse unit-pivot
// pre-elimination (SmithUnitPivotEliminate) before handing what is left to
// its backend. A unit pivot splits off an invariant factor 1 by a pair of
// unimodular operations, and on a sparse matrix the elimination touches
// almost nothing, so the trivial part of the Smith form is stripped at
// negligible cost. The pass is self-limiting (it gives up on a dense or
// unit-free matrix after the first pivot search), hence it is on by
// default for every implementation of Z.
template <typename T> struct use_unit_pivot_preelimination {
  static const bool value = is_implementation_of_Z<T>::value;
};

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
