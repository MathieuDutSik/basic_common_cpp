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
(2026-11, commit 112140f8 there). Over `fmpz_class` / `fmpq_class` the
products written in the natural Eigen way, `G * B.transpose()` or
`A * Ainv * B`, were among the main costs, and the fix was to write them as
`MatrixProduct(MatrixProduct(A, Ainv), B)` and
`MatrixProduct(G, TransposedMat(B))` (`src_matrix/MAT_MatrixFund.h`), which
routes the flint types to `fmpz_mat_mul` / `fmpq_mat_mul`.

That is not an acceptable way to write the code. The expression `a * b * c`
is the one we want to write, and with working expression templates it should
also be the fast one, with no intermediate type to name or to copy.

### Measurements (Apple Silicon, -O3, entries of a few bits)

| product                       | Eigen `A * B.transpose()` | `MatrixProduct(A, TransposedMat(B))` | of which `TransposedMat` |
|-------------------------------|------:|------:|------:|
| fmpz 19x19 * 19x1             |  2.7 us |  1.5 us | 0.06 us |
| fmpz 19x19 * 19x4             | 10.9 us |  2.2 us | 0.20 us |
| fmpz 19x19 * 19x12            | 29.7 us |  3.9 us | 0.45 us |
| fmpz 19x19 * 19x19            | 46.7 us |  5.7 us | 0.87 us |
| fmpq 19x19 * 19x4             | 21.3 us |  8.0 us | 0.27 us |
| fmpq 19x19 * 19x12            | 62.1 us | 12.9 us | 0.76 us |
| fmpq 19x19 * 19x19            |  108 us | 18.7 us | 1.41 us |
| fmpz 5x5 * 5x5                |  1.35 us | 0.37 us | 0.07 us |
| fmpq 5x5 * 5x5                |  2.18 us | 1.23 us | 0.11 us |

The results are equal. So the named product wins by 2x to 8x although it
materializes the transpose, whose copy is under 10% of its cost.

### Why the expression templates do not help

- Eigen's expression templates work at the level of the matrix
  expression: `B.transpose()` is a lazy view and `A * B * C` is
  evaluated as `(A * B) * C` into temporaries (a product is never fused
  lazily). That part is fine. The cost is in the evaluation itself: for a
  scalar type it does not know, Eigen's generic GEBP kernel does one scalar
  `operator*` / `operator+=` call per term. For `fmpq_class` each of these
  also canonicalizes (a gcd). flint's `fmpz_mat_mul` packs small entries
  into machine words, and `fmpq_mat_mul` clears the denominators once and
  multiplies integer matrices. Eigen has no way to know that these kernels
  exist.
- The scalar-level expression templates (gmpxx `__gmp_expr`, the
  `fmpq_product` of `NumberTheoryFlint.h`) do not reach inside Eigen's
  kernel either. Eigen forms `acc += a * b` on its own temporaries, which
  is the unfused FMA described in the section above.

### What a fix has to provide

- `a * b`, `a * b * c`, `a * b.transpose()`, matrix-vector products and,
  ideally, `.dot()` evaluated by the native kernel for the flint types, with
  no change at the call site. The natural route is to specialize Eigen's
  product evaluator (`Eigen::internal::generic_product_impl<..., GemmProduct>`
  and the `GemvProduct` case) for the `fmpz_class` / `fmpq_class` scalars,
  as Eigen does for its BLAS backends. Check whether `fmpz_class` has the
  layout of one `fmpz` (one `slong`); if so an `fmpz_mat_t` window can be
  built on the Eigen storage without copying it, transposes and blocks
  included (strides).
- gmp and flint have different designs, so they need different answers:
  - flint has the matrix kernels; the work is the bridge.
  - gmp has no matrix kernel. For `mpq_class` the same idea as
    `fmpq_mat_mul_cleared` can be written once: scale the rows of A and the
    columns of B to integers, multiply the `mpz_class` matrices, divide once
    per entry. For `mpz_class` the per-term temporary of `acc += a * b` (see
    above) is the first thing to remove, through a specialized kernel using
    `mpz_addmul`.
  - The boost types are a third case (`mpz_int` fuses `acc += a*b` itself).
- It must stay exact and dimension independent. The sizes that matter are
  small (5 to 20, with a few rows of up to a few hundred), so the dispatch
  overhead has to be smaller than these kernels.

### Call sites to revert to `a * b * c` once this works

`MatrixProduct` is now used for this reason in polyhedral_common at
`src_latt/InvariantVectorFamily.h` (projection step of
`inner_span_iterated_shortest_rec`), `src_latt/LatticeRootDecomposition.h`
(`WeylChamberCharacteristicFamily`), `src_latt/LatticeRootSystem.h`
(Dynkin graph), `src_latt/LatticeAutomorphismVinberg.h` (conjugation of the
generators) and `src_group/PolytopeEquiStab.h`
(`T_TranslateToMatrix_ListMat_SHV`). `grep -rn 'MatrixProduct(' src_*` lists
all of them, the older ones included.

### Action items

- [ ] Specialize the Eigen product evaluator for `fmpz_class` / `fmpq_class`
      so that `a * b * c` reaches `fmpz_mat_mul` / `fmpq_mat_mul`, without
      copying when the layout allows it; benchmark on the table above.
- [ ] A product kernel for `mpq_class` (cleared denominators) and for
      `mpz_class` (`mpz_addmul`, no per-term temporary).
- [ ] The same for matrix-vector products and dot products.
- [ ] Then replace the `MatrixProduct(MatrixProduct(a, b), c)` call sites by
      `a * b * c` and remove `MatrixProduct`, or keep it only as an alias.
