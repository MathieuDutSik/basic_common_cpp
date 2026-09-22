// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
//
// Test of underlying_z_ring and underlying_q_field, the rational integers Z
// and the rational numbers Q sitting inside a type, against underlying_ring,
// which answers a different question -- a ring over which the computation can
// run without denominators.
//
// The two agree on every type whose scalars are already rational, which is why
// the distinction stayed invisible for so long. They part ways exactly on the
// algebraic types: for QuadField<mpq_class, d> underlying_ring is Z[sqrt(d)]
// and underlying_z_ring is mpz_class, and for a real algebraic field they are
// the order Z[x] and mpz_class.
//
// What is checked:
//   --- the mapping of every type that has one,
//   --- the absences: underlying_q_field on a ring, and both traits on a type
//       with no rational scalars at all,
//   --- the laws tying the three traits together, which is what keeps a new
//       specialization from drifting away from the others,
//   --- that a value really does travel from a type to its rational integers
//       and back.
//
// Nearly all of it is static_assert, so a failure is a compile error and the
// run only reports that the checks that need values also passed.

// clang-format off
#include "NumberTheoryBoostCppInt.h"
#include "NumberTheoryBoostGmpInt.h"
#include "NumberTheoryCommon.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheory.h"
#include "NumberTheoryQuadField.h"
#include "NumberTheoryRealField.h"
#include "NumberTheoryThreshold.h"
#include "Fp.h"
// clang-format on
#include <string>

// The absence of a trait is a property the code branches on, so it has to be
// detectable rather than a hard error. A specialization whose member typedef
// is invalid to instantiate would not be: it is the specialization itself that
// must not apply, which is what the requires clause on the QuadField case of
// underlying_q_field is for.
template <typename T>
concept HasUnderlyingZRing = requires {
  typename underlying_z_ring<T>::ring_type;
};

template <typename T>
concept HasUnderlyingQField = requires {
  typename underlying_q_field<T>::field_type;
};

template <typename T> using Zof = typename underlying_z_ring<T>::ring_type;
template <typename T> using Qof = typename underlying_q_field<T>::field_type;

// The laws. A field has both traits, its Z is an implementation of Z, its Q is
// an implementation of Q, the Z of its Q is its own Z, and the field of
// fractions of its Z is its Q. Together they pin a new specialization down
// against the ones already there.
template <typename T> constexpr bool CheckField() {
  static_assert(HasUnderlyingZRing<T>, "a field has rational integers");
  static_assert(HasUnderlyingQField<T>, "a field contains Q");
  static_assert(is_implementation_of_Z<Zof<T>>::value,
                "underlying_z_ring lands on an implementation of Z");
  static_assert(is_implementation_of_Q<Qof<T>>::value,
                "underlying_q_field lands on an implementation of Q");
  static_assert(std::is_same_v<Zof<Qof<T>>, Zof<T>>,
                "the rational integers of Q inside T are those of T");
  // The way back up is overlying_field, not underlying_q_field: the rational
  // integers are a ring and so have no Q of their own.
  static_assert(!HasUnderlyingQField<Zof<T>>,
                "the rational integers of T are a ring");
  static_assert(
      std::is_same_v<typename overlying_field<Zof<T>>::field_type, Qof<T>>,
      "the field of fractions of the rational integers of T is Q inside T");
  return true;
}

// A ring has rational integers and no Q, and is its own rational integers when
// it is an implementation of Z.
template <typename T> constexpr bool CheckZRing() {
  static_assert(HasUnderlyingZRing<T>, "a Z ring has rational integers");
  static_assert(!HasUnderlyingQField<T>, "a ring contains no Q");
  static_assert(is_implementation_of_Z<T>::value, "the type is a Z ring");
  static_assert(std::is_same_v<Zof<T>, T>, "a Z ring is its own Z");
  return true;
}

// ---- the plain rational types: the two traits agree with underlying_ring ----

static_assert(CheckZRing<mpz_class>());
static_assert(CheckField<mpq_class>());
static_assert(std::is_same_v<Zof<mpq_class>, mpz_class>);
static_assert(std::is_same_v<Qof<mpq_class>, mpq_class>);
static_assert(std::is_same_v<Zof<mpq_class>,
                             underlying_ring<mpq_class>::ring_type>,
              "on a rational type the two rings coincide");

using Tmpz_b = boost::multiprecision::mpz_int;
using Tmpq_b = boost::multiprecision::mpq_rational;
static_assert(CheckZRing<Tmpz_b>());
static_assert(CheckField<Tmpq_b>());
static_assert(std::is_same_v<Zof<Tmpq_b>, Tmpz_b>);
static_assert(std::is_same_v<Zof<Tmpq_b>, underlying_ring<Tmpq_b>::ring_type>);

using Tcpp_z = boost::multiprecision::cpp_int;
using Tcpp_q = boost::multiprecision::cpp_rational;
static_assert(CheckZRing<Tcpp_z>());
static_assert(CheckField<Tcpp_q>());
static_assert(std::is_same_v<Zof<Tcpp_q>, Tcpp_z>);
static_assert(std::is_same_v<Zof<Tcpp_q>, underlying_ring<Tcpp_q>::ring_type>);

static_assert(CheckZRing<SafeInt64>());
static_assert(CheckField<Rational<SafeInt64>>());
static_assert(std::is_same_v<Zof<Rational<SafeInt64>>, SafeInt64>);
static_assert(std::is_same_v<Qof<Rational<SafeInt64>>, Rational<SafeInt64>>);

static_assert(CheckZRing<int64_t>());
static_assert(CheckField<Rational<int64_t>>());
static_assert(std::is_same_v<Zof<Rational<int64_t>>, int64_t>);
static_assert(std::is_same_v<Zof<int32_t>, int32_t>);
static_assert(std::is_same_v<Zof<int16_t>, int16_t>);
static_assert(std::is_same_v<Zof<int8_t>, int8_t>);

#ifdef ENABLE_FLINT_SUPPORT
static_assert(CheckZRing<fmpz_class>());
static_assert(CheckField<fmpq_class>());
static_assert(std::is_same_v<Zof<fmpq_class>, fmpz_class>);
static_assert(std::is_same_v<Qof<fmpq_class>, fmpq_class>);
// The backends do not cross: a flint type stays on flint.
static_assert(!std::is_same_v<Zof<fmpq_class>, mpz_class>);
#endif

// ---- the quadratic fields: here the two rings differ ----

template <typename Trat, int d> constexpr bool CheckQuad() {
  using T = QuadField<Trat, d>;
  using Tring = QuadField<typename underlying_ring<Trat>::ring_type, d>;
  static_assert(CheckField<T>());
  // underlying_ring keeps sqrt(d), underlying_z_ring leaves it behind.
  static_assert(std::is_same_v<typename underlying_ring<T>::ring_type, Tring>,
                "the fraction-free ring of Q(sqrt(d)) is Z[sqrt(d)]");
  static_assert(std::is_same_v<Zof<T>, Zof<Trat>>,
                "the rational integers of Q(sqrt(d)) are those of the base");
  static_assert(!std::is_same_v<Zof<T>, typename underlying_ring<T>::ring_type>,
                "the two rings of a quadratic field are different types");
  static_assert(std::is_same_v<Qof<T>, Qof<Trat>>,
                "the rationals of Q(sqrt(d)) are those of the base");
  // Z[sqrt(d)] is a ring: rational integers, but no Q.
  static_assert(HasUnderlyingZRing<Tring>, "Z[sqrt(d)] has rational integers");
  static_assert(!HasUnderlyingQField<Tring>, "Z[sqrt(d)] contains no Q");
  static_assert(std::is_same_v<Zof<Tring>, Zof<Trat>>);
  return true;
}

static_assert(CheckQuad<mpq_class, 2>());
static_assert(CheckQuad<mpq_class, 3>());
static_assert(CheckQuad<mpq_class, 5>());
static_assert(CheckQuad<Tmpq_b, 5>());
static_assert(CheckQuad<Rational<SafeInt64>, 5>());
static_assert(std::is_same_v<Zof<QuadField<mpq_class, 5>>, mpz_class>);
static_assert(std::is_same_v<Qof<QuadField<mpq_class, 5>>, mpq_class>);
// The base backend is followed, so a boost base does not fall back on gmp.
static_assert(std::is_same_v<Zof<QuadField<Tmpq_b, 5>>, Tmpz_b>);
#ifdef ENABLE_FLINT_SUPPORT
static_assert(std::is_same_v<Zof<QuadField<fmpq_class, 5>>, fmpz_class>);
#endif

// ---- the real algebraic field: the coefficients are fixed by the build ----

static_assert(CheckField<RealField<1>>());
static_assert(std::is_same_v<Zof<RealField<1>>, Tint_real_field>);
static_assert(std::is_same_v<Qof<RealField<1>>, Trat_real_field>);
static_assert(
    std::is_same_v<underlying_ring<RealField<1>>::ring_type, RealRing<1>>,
    "the fraction-free ring of the field is the order Z[x]");
static_assert(!std::is_same_v<Zof<RealField<1>>, RealRing<1>>,
              "the two rings of a real algebraic field are different types");
// Z[x] is a ring: rational integers, but no Q.
static_assert(HasUnderlyingZRing<RealRing<1>>);
static_assert(!HasUnderlyingQField<RealRing<1>>);
static_assert(std::is_same_v<Zof<RealRing<1>>, Tint_real_field>);
// Distinct field indices answer the same, the coefficients not depending on
// which field was registered.
static_assert(std::is_same_v<Zof<RealField<2>>, Zof<RealField<1>>>);

// ---- the types with no rational scalars: both traits absent ----

static_assert(!HasUnderlyingZRing<double>, "double has no rational integers");
static_assert(!HasUnderlyingQField<double>);
static_assert(!HasUnderlyingZRing<float>);
static_assert(!HasUnderlyingQField<float>);
// A finite field: it is neither an implementation of Z nor one of Q, and no
// copy of either sits inside it.
static_assert(!HasUnderlyingZRing<Fp<int64_t, 2147389441>>);
static_assert(!HasUnderlyingQField<Fp<int64_t, 2147389441>>);
// Inexact, so neither Z nor Q is represented exactly.
static_assert(!HasUnderlyingZRing<ThresholdField<1>>);
static_assert(!HasUnderlyingQField<ThresholdField<1>>);

// ---- the values follow the types ----

static int n_error = 0;

static void check(bool test, std::string const &name) {
  if (test) {
    std::cerr << "PASS: " << name << "\n";
  } else {
    std::cerr << "FAIL: " << name << "\n";
    n_error++;
  }
}

// An integer of T travels to the rational integers of T and back unchanged,
// and a value that is not a rational integer is refused on the way down.
template <typename T>
static void RoundTrip(std::string const &name, T const &integral,
                      T const &non_integral, bool has_non_integral) {
  using Tz = Zof<T>;
  Tz down = UniversalScalarConversion<Tz, T>(integral);
  T up = UniversalScalarConversion<T, Tz>(down);
  check(up == integral, name + ": an integer travels down to Z and back");
  check(down == UniversalScalarConversion<Tz, int>(7),
        name + ": and it is the same integer");
  if (has_non_integral) {
    check(!UniversalScalarConversionCheck<Tz, T>(non_integral).has_value(),
          name + ": a value outside Z is refused on the way down");
  }
}

int main() {
  try {
    using Tq2 = QuadField<mpq_class, 2>;
    using Tq5 = QuadField<mpq_class, 5>;
    RoundTrip<mpq_class>("mpq_class", mpq_class(7), mpq_class(3, 4), true);
    RoundTrip<Tmpq_b>("mpq_rational", Tmpq_b(7), Tmpq_b(3) / Tmpq_b(4), true);
    RoundTrip<Rational<SafeInt64>>("Rational<SafeInt64>",
                                   Rational<SafeInt64>(7),
                                   Rational<SafeInt64>(3) /
                                       Rational<SafeInt64>(4),
                                   true);
    RoundTrip<Tq2>("Qsqrt2", Tq2(mpq_class(7), mpq_class(0)),
                   Tq2(mpq_class(0), mpq_class(1)), true);
    RoundTrip<Tq5>("Qsqrt5", Tq5(mpq_class(7), mpq_class(0)),
                   Tq5(mpq_class(0), mpq_class(1)), true);
#ifdef ENABLE_FLINT_SUPPORT
    RoundTrip<fmpq_class>("fmpq_class", fmpq_class(7),
                          fmpq_class(3) / fmpq_class(4), true);
#endif
    if (n_error > 0) {
      std::cerr << "Test_UnderlyingTraits: " << n_error << " error(s)\n";
      throw TerminalException{1};
    }
    std::cerr << "Normal termination of Test_UnderlyingTraits\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_UnderlyingTraits\n";
    exit(e.eVal);
  }
}
