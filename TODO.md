# TODO / performance notes

## Scalar FMA `acc += a * b`: intermediate value vs. direct

Investigation (2026-07) into the cost of the fused-multiply-add pattern
`acc += a * b` that dominates the inner loops of the matrix kernels
(`DeterminantMat*`, `TMat_Inverse_destroy`, ...) and the jet arithmetic
(`jet_number.h`). Three code paths were benchmarked per type:

- `fused`   : `acc += a * b;`
- `scratch` : `prod = a * b; acc += prod;`   (one hoisted `prod`, reused across
              iterations; a GMP-backed temporary reuses its limb buffer on
              reassignment instead of allocating/freeing every iteration)
- `inplace` : `prod = a; prod *= b; acc += prod;`  (never constructs the fresh
              value that binary `operator*` returns)

64-term dot products, small operands (the regime of the jet coefficients),
`-O2`. Times in ms, lower is better; the winner is in **bold**.

### Integer types
| type                 | fused | scratch | inplace |
|----------------------|------:|--------:|--------:|
| `mpz_class` (gmpxx)  |   240 |  **77** |      97 |
| boost `mpz_int`      |**75** |      83 |     104 |
| boost `cpp_int`      |    98 |  **55** |      83 |

### Rational types
| type                     | fused | scratch | inplace |
|--------------------------|------:|--------:|--------:|
| `mpq_class` (gmpxx)      |  1672 |**1027** |    1060 |
| boost `mpq_rational`     |  1671 |**1022** |    1057 |
| boost `cpp_rational`     |  1580 |    1565 |    1585 | (neutral)

### Composed types
| type                        | fused | scratch | inplace |
|-----------------------------|------:|--------:|--------:|
| `Rational<mpz_class>`       |  1893 |    1889 |**1820** |
| `Rational<long>`            |   173 |     173 |     173 | (neutral)
| `QuadField<mpq_class, 3>`   |   866 |     985 | **811** |

### Findings

1. **gmpxx does NOT fuse `acc += a*b`.** For `mpz_class`/`mpq_class` it
   materializes a freshly-allocated temporary for `a*b` every iteration, so the
   `scratch` form (reusing one `prod`) is a large win: 3.1x for `mpz_class`,
   1.6x for `mpq_class`. `prod = a*b` is fast because `a*b` is a lazy
   `__gmp_expr` evaluated in place (`mpz_mul(prod, a, b)`).

2. **boost `mpz_int` is the opposite.** Its expression templates DO fuse
   `acc += a*b` (fused is the fastest path for it); `scratch`/`inplace` make it
   slower. So `mpz_class` and `mpz_int` want OPPOSITE code paths -- there is no
   single best FMA form across types. `cpp_int` likes `scratch`;
   `cpp_rational` is neutral.

3. **Value-returning composite types (`Rational`, `QuadField`) want `inplace`,
   not `scratch`.** Their binary `operator*` builds a fresh value
   (`Rational z; ...; return z;`), so an outer `prod` cannot help (`scratch` is
   neutral for `Rational`, and 13% WORSE for `QuadField`). The compound `*=`
   avoids the fresh value and is ~4% (`Rational<mpz>`) / ~7% (`QuadField`)
   faster than fused. The ceiling is modest because the inner mpz/mpq
   temporaries and GCD reductions dominate, not the wrapper allocation.

### What is actually used in this repo

The hot type is `mpq_class` (jet coefficients: `jet<mpq_class, N>`). For it
`scratch` is unambiguously best. That win is already applied inside
`jet_number.h` `operator*` / `inverse` (their inner loops multiply bare
`mpq_class` coefficients). This is safe because those loops only ever see
`mpq_class`.

### Action items

- [ ] A `scalar_fma(acc, a, b)` helper gated by a trait (e.g.
      `fma_prefers_scratch<T>`, default false; true for `mpz_class`,
      `mpq_class`, boost `cpp_int`/`mpq_rational`/`cpp_rational` where scratch
      helps; explicitly FALSE for boost `mpz_int`, which fuses natively). This
      is the only clean way to serve both `mpz_class` (scratch) and `mpz_int`
      (fused) from one generic call site. Keep it minimal -- do not build a
      general expression-template framework.
- [ ] Do NOT apply `scratch` blindly in the generic matrix kernels: there
      `T` may be a jet / `QuadField` / native, for which scratch is
      neutral-to-harmful. If the kernels are optimized, use `inplace` (`*=`) or
      gate on the trait.
- [ ] (Lower priority) `rational.h` / `NumberTheoryQuadField.h`: binary
      `operator*` constructs a fresh value; an expression-template layer (lazy
      `MulExpr` evaluated in `operator=`) would let `prod = a*b` reuse the
      destination's buffers, matching `mpz_class`. Possible but modest payoff
      (~4-5%) for these types -- their cost is dominated by inner arithmetic /
      GCDs, not the wrapper allocation. `Rational<mpz_class>` is ~12x slower
      than `mpq_class` regardless; prefer `mpq_class` where possible.

## Matrix products of multiprecision types: `a * b * c` should be the fast path

Found while speeding up the Weyl-chamber canonical form of polyhedral_common
(2026-11, commit 112140f8 there): over `fmpz_class` / `fmpq_class` the
products written the natural Eigen way, `G * B.transpose()` or
`A * Ainv * B`, were among the main costs, and they had been rewritten as
`MatrixProduct(MatrixProduct(A, Ainv), B)` to reach `fmpz_mat_mul` /
`fmpq_mat_mul`.

### What was going on

The gap was never about the transpose or about the intermediate of
`A * B * C`: Eigen's `B.transpose()` is a free view, and `A * B * C` is
evaluated as `(A * B) * C` into a temporary whatever the type, exactly as
the nested `MatrixProduct`. It was the kernel: `A * B` ran Eigen's generic
GEMM, one scalar `acc += a * b` per term (a call, the small/large check of
the fmpz, a gcd per term for fmpq), while `MatrixProduct` called flint,
which checks the bit sizes once and runs the inner loop on machine words
(fmpz), or clears the denominators once and does one integer product
(fmpq). For `int`, `mpz_class`, `mpq_class` both were the same code.

### Done (2026-11): `a * b` itself goes to flint

`src_matrix/EigenFlintProduct.h`, included by `MatrixTypes.h` with
`ENABLE_FLINT_SUPPORT`: partial specializations, constrained to
`is_flint_type<T>`, of `Eigen::internal::general_matrix_matrix_product`
and `general_matrix_vector_product`, the class templates through which
Eigen evaluates every matrix-matrix and matrix-vector product (as its BLAS
backend does for float/double). They view the Eigen storage as a flint
matrix without copying (`fmpz_class` is one `fmpz`, a flint matrix is
row-major entries with a row stride, so a column-major Eigen matrix is
viewed as its transpose and the product computed as `C^T = B^T A^T`); an
operand in the other order is copied transposed. `EigenProductThreshold.h`
makes Eigen's small-product cutoff (`EIGEN_GEMM_TO_COEFFBASED_THRESHOLD`)
depend on the scalar so that the small products go there too.
`MatrixProduct` and `FlintProductMatrix` are removed, and the call sites in
polyhedral_common are back to `a * b * c`. Tested by
`src_matrix/Test_EigenFlintProduct`.

Measurements (Apple Silicon, -O3, entries in [-7, 7], us):

| product        | Eigen generic | old `MatrixProduct` | now `*` |
|----------------|------:|------:|------:|
| fmpz 5x5 * 5x5     |  1.20 |  0.32 | 0.20 |
| fmpz 19x19 * 19x19 |  55.6 |  5.09 | 3.87 |
| fmpq 19x19 * 19x19 |  96.3 |  12.0 | 10.3 |
| fmpz `A * B * C` 19 | 112  |  10.2 | 7.7 |
| fmpz 60x60 * vector |  30.3 |  12.2 (n x 1) | 4.7 |
| fmpq 60x60 * vector |  45.6 |  31.9 (n x 1) | 20.0 |

### Remaining

- [ ] With Eigen 3.3 (the CI version) the cutoff is the literal 20, not the
      macro: the products with rows + cols + depth < 20 (about 6x6 and
      below) stay on the generic coefficient loop there. Correct, slower.
- [ ] `.dot()` and the 1x1 product `v.transpose() * w` go through Eigen's
      coefficient reduction, not through these kernels.
- [ ] A product kernel for `mpq_class` (cleared denominators) and for
      `mpz_class` (`mpz_addmul`, no per-term temporary): gmp has no matrix
      kernel, so this has to be written. `mpz_class` 19x19 takes 288 us
      in the generic kernel against 3.9 us for `fmpz_class` above.
- [ ] Column-major vs row-major storage of `MyMatrix`: to be discussed.
      Not needed for the products; it matters for the flint routines that
      work in place (Hermite normal form), which see the transpose.
