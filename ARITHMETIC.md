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

### FLINT types (optional, `ENABLE_FLINT_SUPPORT`)

  * **fmpz_class** -- wrapper around the `fmpz_t` integer of the
    [flint](https://flintlib.org) library. Header:
    `src_number/NumberTheoryFlint.h`.
  * **fmpq_class** -- wrapper around `fmpq_t`, the flint rational. Same
    header.

The wrappers expose the interface surface of `mpz_class` / `mpq_class`
(operators, `get_num` / `get_den`, iostream, the trait specializations,
`TYPE_CONVERSION` overloads including conversions to and from the gmp
types), so the generic matrix code accepts them unchanged. The point of
flint over gmp is the small-integer optimization: an `fmpz` is a single
word holding either the value itself (below 62 bits) or a pointer to an
mpz, so matrices of typical small entries never touch the allocator. The
multiply-accumulate goes through the native fused `fmpz_addmul` /
`fmpq_addmul` calls (see `is_fma_prefered`).

The support is compiled only when `ENABLE_FLINT_SUPPORT` is defined; in
`src_matrix` this is done with `make ENABLE_FLINT_SUPPORT=1` (the
`FLINT_INCDIR` / `FLINT_LINK` variables override the include directory and
the link flags, the defaults assume flint installed next to gmp). The
programs then accept the arithmetic names `flint_integer` and
`flint_rational` next to `integer` / `rational`, and
`src_matrix/Bench_matrix_arithmetic` benchmarks the two families on
identical inputs.

Three mechanisms stack on top of the plain wrappers:

  * **Deferred products.** `a * b` between two flint values returns a
    lightweight proxy, and the consuming operation picks the fused call:
    `acc += a*b` becomes `fmpz_addmul`, `x = a*b - c*d` a mul/submul
    chain, with no temporaries. Any other use converts the proxy back to
    the concrete class, so generic code and Eigen compile unchanged. Same
    caveat as the gmpxx expression templates: do not store `auto p = a*b;`
    beyond the full expression.
  * **Native matrix backend** (`src_matrix/MAT_MatrixFlint.h`). The named
    generic operations dispatch `MyMatrix<fmpz_class>` /
    `MyMatrix<fmpq_class>` to `fmpz_mat` / `fmpq_mat`, which carry the
    machine-word and multimodular algorithms: `MatrixProduct`,
    `DeterminantMat`, the row and column Hermite normal forms (all
    variants; the column form goes through the transpose, whose convention
    is the exact mirror), `SmithNormalFormInvariant`, `Inverse` (both
    types; the inverse is unique so the result is identical to the generic
    one) and `NullspaceTrMat` / `NullspaceMat` (through `fmpq_mat_rref`;
    the reduced row echelon form is unique and the kernel extraction uses
    the same formula, so the basis is identical to the generic one). The
    call sites stay fully generic; the marshalling is O(n^2) word copies,
    measured below 3% of the cheapest routed operation. Measured against
    the gmp types on identical inputs: small-entry 120x120 products
    73ms -> 0.6ms (integer) and 273ms -> 0.9ms (rational), 60x60
    determinants 5.8ms -> 0.45ms and 18.6ms -> 0.46ms, 32x32
    HNF-with-transform 1.13s -> 0.66ms (the modular HNF has no coefficient
    explosion), 50x50 rational inverse 62ms -> 5ms and 60x90 rational
    nullspace 54ms -> 2ms.
  * **Modulo-D Hermite normal form** (generic, in `MAT_MatrixInt.h`).
    `HermiteNormalFormModD_or_none` implements the Domich-Kannan-Trotter
    scheme for ANY exact euclidean domain opting in through
    `use_hnf_mod_D<T>` (currently `mpz_class` and `fmpz_class`; a future
    Z[i] can opt in): all intermediate entries stay bounded by the
    determinant D of a row selection, and a final back-substitution
    certificate proves the result exact (falling back to the generic
    kernel otherwise). All the Hermite entry points dispatch to it: the
    H-only forms directly, the column forms through the transpose, and
    the (U, H) pair forms for a square nonsingular M, where the unique
    transform U = H M^{-1} is recovered by one inversion over the
    overlying field. For mpz_class this took the 32x32 HNF from ~1.1s to
    ~4ms and makes 100x100 (formerly out of reach) run in ~0.3s.

  * **Sparse splitting-pivot pre-elimination for the Smith form**
    (generic, in `MAT_MatrixInt.h`). An entry u = M(p, q) that divides
    every entry of its own row and of its own column is a *splitting
    pivot*: both eliminations are exact, M becomes diag(u) + core, and the
    cokernel is a direct sum. `SmithUnitPivotEliminate` peels those off
    before the backend sees anything. A unit is the special case u = +-1,
    which needs no verification and covers a torsion-free matrix; a
    general splitting pivot is what carries the torsion, and it is
    recognized by a divisibility test with nothing to choose or tune.
    Because the cokernel splits, the elementary divisors of the whole are
    the union of those of the parts, and the divisibility chain is
    restored at the end by the gcd / lcm cascade -- the Smith form of the
    diagonal matrix carrying all the factors. With only unit pivots that
    cascade is a no-op and the assembly is a concatenation.

    The pass works on a sparse image of the matrix and chooses its pivot
    by the Markowitz criterion, minimizing (|row p| - 1) * (|col q| - 1);
    the search is symmetric in rows and columns, which matters because a
    boundary matrix of a graph has a dozen entries per row and exactly two
    per column, so its cheap pivots are only visible from the column side.
    It runs the search for units first and falls back to the general test
    only when that finds nothing -- an ordering by cost, not a preference
    between pivots. `use_unit_pivot_preelimination<T>` gates the whole
    thing, on by default for every implementation of Z, and
    `IsUnitForPivot` is the customization point for a ring with more units
    than +-1 (Z[i]).

    The pass stops on two conditions: no splitting pivot is left, or the
    fill has made what remains dense, which is the point where handing
    over to a dense backend is the better move. That density is derived
    rather than tuned: a sparse entry costs a row index, a column index
    and a value against the single value of the dense form, so the sparse
    representation stops paying for itself around one third. A bound on
    the fill of a single pivot used to sit alongside it and was removed
    after being measured to be inert -- 2^20 and 2^40 gave bit-identical
    results, the elimination being stopped by the exhaustion of its
    pivots, never by the fill of one of them. The one quantity left is
    how many rows and columns the Markowitz search looks at, and it is a
    search effort where more is not better: an unbounded window is an
    exact Markowitz minimum, which on a 11568 x 12119 matrix costs 121
    seconds against 2.3 and leaves a worse core, the exact minimum of a
    greedy criterion not being a global optimum.

    A sparse input path avoids the dense form altogether:
    `SmithNormalFormInvariant_sparse` in `MAT_MatrixIntSparse.h` takes a
    `MySparseMatrix`, runs the same elimination on it, and expands only
    the core, if any is left. This matters because the whole dense
    pipeline costs time proportional to the number of CELLS: parsing the
    file, constructing that many number objects, scanning for the density
    and building the sparse image. On an 11568 x 12119 boundary matrix the
    elimination itself is 1.7 s of a 7.9 s dense-input run, so about 78%
    of the work was handling a dense form of a matrix that is not dense.
    `ConvertMatrixDenseToSparse` rewrites a dense matrix file into the
    sparse format, streaming the entries so that a matrix too large to
    hold densely can still be converted, and
    `CI_tests/SNF_computations` holds one chain complex in that form
    (582 MB of dense text becomes 2.9 MB).

    This is the regime where both dense backends are at their worst, their
    cost following the dimensions rather than the number of nonzero
    entries. On the boundary matrices of a chain complex (0.06% to 2%
    dense, entries almost all +-1, invariant factors nearly all 1) the
    computation collapses in the pre-elimination, the backend receiving
    only what is left:

        matrix        backend alone      with the pre-elimination
        1477x126      75 ms mpz          2.4 ms mpz
                      691 ms flint       2.6 ms flint
        870x5253      21.0 s mpz         115 ms mpz / 50 ms flint
        6285x1477     48.9 s mpz         169 ms mpz / 198 ms flint

    The invariant factors are identical to the ones of the backend alone on
    all three.

    The same machinery extends to the Smith invariant factors:
    `SmithNormalFormInvariantModD_or_none` runs the Kannan-Bachem
    alternation of modulo-D Hermite reductions (each pass a certified
    one-sided unimodular equivalence) until the matrix is diagonal, then
    reads the invariants off the gcd/lcm cascade of the diagonal. It is
    gated by its own trait `use_snf_mod_D<T>`, which currently NO type
    sets: unlike the row-only Hermite elimination, the generic Smith
    kernel reduces from both sides, shows no coefficient explosion, and
    beats the alternation by ~1.2-2x at every size probed (16..96,
    spreads 10..10^6). The machinery is kept for a ring whose generic
    Smith kernel does blow up; `fmpz_class` uses the native
    `fmpz_mat_snf` instead.

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
