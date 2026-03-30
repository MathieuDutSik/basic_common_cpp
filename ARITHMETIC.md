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
  * `underlying_ring<T>` -- the ring of integers inside a field type.

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
