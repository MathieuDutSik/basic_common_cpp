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
