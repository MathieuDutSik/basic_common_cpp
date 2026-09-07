Arithmetics
===========

Many of the programs in this repository allow different arithmetic types to be
used. This is both for functional and speed reasons. The types are always
template parameters, which allows compile-time selection and zero-overhead
abstraction.

## Type system overview

The code relies on a traits-based system defined in `src_number/TemplateTraits.h`
to classify types at compile time. The key traits are:

  * `is_euclidean_domain<T>` -- whether GCD and Euclidean division are available.
  * `is_implementation_of_Z<T>` -- whether the type represents the integers Z.
  * `is_implementation_of_Q<T>` -- whether the type represents the rationals Q.
  * `overlying_field<T>` -- the field of fractions of an integer type.
  * `underlying_ring<T>` -- a ring inside a field type over which the
    computation can be run without denominators. It is not canonical and it
    need not be the ring of integers: for a real algebraic field it is the
    order spanned by the powers of the generator (see `RealRing` below).

These traits drive `static_assert` checks and `if constexpr` / SFINAE dispatch
throughout the library.

## Type conversion

Conversions between arithmetic types go through the function
`UniversalScalarConversion<Tout, Tin>(val)` and the lower-level overloads of
`TYPE_CONVERSION(stc<Tin>, Tout&)` defined in `src_number/TypeConversion.h`.
This mechanism avoids implicit narrowing and provides a single point of
control for all inter-type conversions.

## Functional types (algebraic number fields)

### Quadratic fields -- `QuadField<T, d>`

Defined in `src_number/NumberTheoryQuadField.h`.
An element of **Q(sqrt(d))** is stored as a pair `(a, b)` representing
`a + b * sqrt(d)`, where `T` is the coefficient type (typically `mpq_class`).

Pre-defined shorthands used in the command-line interface:

  * **Qsqrt2** -- `QuadField<mpq_class, 2>`, the field **Q(sqrt(2))**.
  * **Qsqrt5** -- `QuadField<mpq_class, 5>`, the field **Q(sqrt(5))**.

To add a new quadratic field such as **Q(sqrt(6))**, add a new branch in
`src_number/NumberTheoryMultiType.h` following the existing pattern.

Elements are written in the input format `a+b*x`, for instance `1`, `1+x`,
`x/4`, `-3-3x/4`, `3+(3*4)*x^3`. No spaces are allowed inside an entry.

#### The underlying ring

`underlying_ring<QuadField<T, d>>::ring_type` is `QuadField<Tring, d>` with
`Tring` the underlying ring of the base, so `QuadField<mpq_class, d>` has
`QuadField<mpz_class, d>` -- that is `Z[sqrt(d)]` -- as its ring. No separate
class is needed: the arithmetic, the ordering, the input/output, the hashing
and the serialization are the same over either base, and the traits already
follow the base type, `is_ring_field<QuadField<T, d>>` being
`is_ring_field<T>`.

`Z[sqrt(d)]` is not the ring of integers of the field: for `d = 1 mod 4` that
is the strictly larger `Z[(1+sqrt(d))/2]`, which the `(a, b)` layout over
`(1, sqrt(d))` cannot represent. It is a ring for every `d` since `sqrt(d)` is
a root of the monic `X^2 - d`, so no monicity check arises here, unlike
`RealRing`.

What changes with the base is the division. Over either base it is
`x / y = x * conj(y) / N(y)` with the norm `N(y) = c^2 - d e^2` a scalar of the
base, and the two coordinates of `x * conj(y)` are divided by it. Over a field
those divisions are exact by construction. Over a ring the quotient lies in
the ring exactly when the norm divides both coordinates, and **when it does
not the division emits an error and throws `TerminalException{1}`** rather
than truncate. It used to truncate silently: `(1+sqrt(5))/2` came back as `0`.

The same `conj(s)/N(s)` gives the canonicalization of a vector inside the ring
(`ScalarCanonicalizationVectorRing`), which is what keeps the coefficients
from growing where the ring has no gcd to reduce a content with.

Measured on cyclic polytopes over `Q(sqrt(5))`, running the dual description
over the ring instead of the field: about 2.3x for beneath-and-beyond, 2.5x
for lrs, and 806x for cdd on the larger case, the double description being the
place where unbounded ray growth was the whole cost.

### General real algebraic fields -- `RealField<i_field>`

Defined in `src_number/NumberTheoryRealField.h`.
For algebraic numbers of degree > 2, the library uses a more general
representation. An element of a number field of degree `d` is stored as a
`std::vector<T>` of `d` rational coefficients in the power basis
`{1, alpha, alpha^2, ..., alpha^(d-1)}`.

Arithmetic operations (especially division) require solving linear systems, so
this type is slower than `QuadField`.

The field is specified at runtime via a descriptor file containing:
  1. The degree `d`.
  2. The `d+1` coefficients of the minimal polynomial.
  3. A `double` approximation of the real root.
  4. A list of rational lower/upper bound pairs used for sign determination
     via continued-fraction approximants.

Usage on the command line: **RealAlgebraic=FileDesc**, where `FileDesc` is the
path to the descriptor file.

A `HelperClassRealField<T>` object is constructed from the file and stored in
a global registry (`list_helper`), keyed by a compile-time integer index
`i_field`. The `RealField<i_field>` class then looks up its helper at
construction time.

### The underlying ring -- `RealRing<i_field>`

`underlying_ring<RealField<i_field>>::ring_type` is `RealRing<i_field>`, the
order `Z[x] = sum_{0 <= i < d} Z x^i` spanned by the powers of the generator.
It is not the ring of integers of the field: no integral closure is computed,
and another generator gives another ring.

`Z[x]` is a ring exactly when the minimal polynomial of `x` is monic, that is
when the coefficient of `x^d` is 1 and `x` is an algebraic integer. **Using
`RealRing` over a descriptor file whose minimal polynomial is not monic emits
an error and throws `TerminalException{1}`.** A non-monic description is easy
to repair: replacing the generator `x` by `c*x` for a suitable integer `c`
makes the minimal polynomial monic without changing the field. For instance
`sin(2*pi/7)`, of minimal polynomial `64 X^6 - 112 X^4 + 56 X^2 - 7`, becomes
`2*sin(2*pi/7)` of minimal polynomial `X^6 - 7 X^4 + 14 X^2 - 7`.

A `RealRing` element is a `d`-tuple of integers with no denominator, so the
gcd normalization that every `RealField` operation performs disappears:
addition and multiplication are plain integer polynomial operations followed
by the reduction rows. Measured on the cubic field of discriminant 49
(`CI_tests/RealAlgebraicField/CubicFieldDisc_49`) with
`src_number/Bench_real_ring`, the ring is about 4x faster on determinants and
1.3x to 2x faster on matrix products. The sign determination is unchanged,
since it uses the same approximant ladder.

Division is where the ring differs in kind from the field. `a / b` is the
solution of `M(b) v = a` with `M(b)` the integer matrix of the multiplication
by `b`, solved fraction-free by `SolveIntegralFractionFree` (Zhou & Jeffrey,
"Fraction-free matrix factors: new forms for LU and QR factors", 2008), so the
computation stays over `Z`. When the solution is not integral the quotient is
not an element of `Z[x]`, there is nothing to fall back on, and the division
emits an error and throws `TerminalException{1}`.

The generic code that runs a computation over `underlying_ring<T>` therefore
sees, for a real algebraic field, a ring that is neither a field nor a
Euclidean domain. Three dispatches in `src_matrix` account for it:
`ScalarCanonicalizationVector` / `ScalarCanonicalizationMatrix` normalize
through the overlying field and scale back into the ring,
`RemoveFractionMatrixPlusCoeffRing` gains the same second regime its vector
counterpart already had, and `SubsetRankOneSolver` uses the
`SubsetRankOneSolver_RingOverField` variant.

`CI_tests/RealAlgebraicField/run_test.sh` covers all of it: the field, the
ring, the field/ring consistency of every one of those paths
(`src_matrix/Test_RingConsistency`) and the timing comparison of the two
(`src_number/Bench_real_ring`). It is run by the number theory and the matrix
CI workflows.

## Speed types (rational and integer implementations)

### GMP types (recommended default)

  * **mpq_class** -- GMP rational type. The default choice for most programs
    since it is well tested and has the fewest issues. Header:
    `src_number/NumberTheoryGmp.h`.
  * **mpz_class** -- GMP arbitrary-precision integer type. Same header.

### Boost.Multiprecision types

  * **boost::multiprecision::cpp_rational** -- Pure C++ rational type from
    Boost.Multiprecision. Header-only (no GMP dependency), but slower than
    `mpq_class`. Header: `src_number/NumberTheoryBoostCppInt.h`.
  * **boost::multiprecision::cpp_int** -- Same library, integer variant.
  * **boost::multiprecision::mpq_rational** -- Boost wrapper around GMP
    rationals. Header: `src_number/NumberTheoryBoostGmpInt.h`.
  * **boost::multiprecision::mpz_int** -- Boost wrapper around GMP integers.
    Same header.

### Template rational -- `Rational<T>`

Defined in `src_number/rational.h`.
A rational number built from an arbitrary integer type `T`, stored as a
numerator/denominator pair with GCD reduction. Typical instantiation:
`Rational<int64_t>` or `Rational<SafeInt64>`.

### Safe bounded integers -- `SafeInt64`

Defined in `src_number/NumberTheorySafeInt.h`.
A wrapper around `int64_t` that checks for overflow on every arithmetic
operation. If the result would exceed the safe bounds (`MAX_INT64_PROD` for
products, `MAX_INT64_SUM` for sums), an exception is thrown. This guarantees
that no silent overflow produces a wrong result, at the cost of some runtime
overhead.

  * **SafeInt64** -- for integer computations in `int64_t` with overflow
    detection.
  * **Rational\<SafeInt64\>** -- rational arithmetic with the same overflow
    safety on the underlying integer operations.

### Finite prime fields -- `Fp<T, P>`

Defined in `src_number/Fp.h`.
A compile-time prime field **F_p** where `P` is a template parameter. Elements
are stored as a single integer of type `T`, kept reduced modulo `P` after every
operation. Division uses the extended Euclidean algorithm for modular inversion.

## p-adic numbers

Defined in `src_number/NumberTheoryPadic.h`.
The library includes support for computations with p-adic numbers, represented
to a fixed precision (degree `d`) as a vector of digits in `{0, ..., p-1}`
together with a valuation exponent. Operations include addition, multiplication,
and inversion (via iterative lifting or extended GCD).

This is used in particular for deciding local square classes
**Q_p\* / (Q_p\*)\^2**.

## Additional utilities

  * **Continued fractions** (`src_number/fractions.h`) -- Compute continued
    fraction expansions and their convergents for rational numbers.
  * **Integer factorization** (`src_number/factorizations.h`) -- Trial division
    and Pollard's rho algorithm for factoring integers.
  * **Quadratic residues** (`src_number/quadratic_residue.h`) -- Computation of
    Legendre symbols and quadratic residues modulo a prime.
  * **GCD and extended GCD** (`src_number/NumberTheoryGeneric.h`) -- Generic
    Euclidean algorithm (`GenericGcd`) and extended GCD (`ComputePairGcdDot`)
    for any Euclidean domain type.
  * **Quotient and residue** (`src_number/QuoIntFcts.h`,
    `src_number/ResidueQuotient.h`) -- Euclidean division primitives `QuoInt`
    and `ResInt`.

## Runtime arithmetic selection

The function `process_by_numeric_type` in `src_number/NumberTheoryMultiType.h`
dispatches on a string argument to instantiate a templated function with the
appropriate type. Currently recognized values:

  * `"rational"` -- `mpq_class`
  * `"Qsqrt2"` -- `QuadField<mpq_class, 2>`
  * `"Qsqrt5"` -- `QuadField<mpq_class, 5>`
  * `"RealAlgebraic=<file>"` -- `RealField<1>` with the field described in
    `<file>`

The default choice for most applications is **mpq_class**.
