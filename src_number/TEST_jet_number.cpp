// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "MAT_Matrix.h"
#include "jet_number.h"
// clang-format on

// Standalone unit tests for the jet<T, N> numeric type and the type-aware pivot
// selection it needs. Build with
//   make CHOICE_COMPILATION= TEST_jet_number
// and, for the internal invariant checks,
//   make CHOICE_COMPILATION=-DSANITY_CHECK TEST_jet_number

template <typename T, int N> using J = jet<T, N>;

int main() {
  using T = mpq_class;
  int n_error = 0;
  auto check = [&](bool ok, std::string const &name) {
    std::cerr << (ok ? "ok   " : "FAIL ") << name << "\n";
    if (!ok)
      n_error++;
  };

  // Arithmetic identities at order 3.
  {
    constexpr int Nd = 3;
    J<T, Nd> t = J<T, Nd>::var();
    J<T, Nd> a = J<T, Nd>(2) + J<T, Nd>(3) * t; // 2 + 3 t
    J<T, Nd> b = J<T, Nd>(1) - t;               // 1 - t
    J<T, Nd> s = a + b;
    check(s[0] == T(3) && s[1] == T(2) && s[2] == T(0) && s[3] == T(0),
          "sum 2+3t plus 1-t");
    J<T, Nd> p = a * b; // 2 + t - 3 t^2
    check(p[0] == T(2) && p[1] == T(1) && p[2] == T(-3) && p[3] == T(0),
          "product (2+3t)(1-t)");
    J<T, Nd> inv = inverse(b); // 1 + t + t^2 + t^3
    check(inv[0] == T(1) && inv[1] == T(1) && inv[2] == T(1) && inv[3] == T(1),
          "geometric series 1/(1-t)");
    J<T, Nd> q = a / b;
    check(q * b == a, "division round-trip (a/b)*b == a");
    T t0(1);
    t0 /= T(10);
    check((a + b).eval(t0) == a.eval(t0) + b.eval(t0), "eval homomorphism +");
  }

  // Leading-coefficient ordering / positivity.
  {
    constexpr int Nd = 2;
    J<T, Nd> t = J<T, Nd>::var();
    J<T, Nd> zero;
    check(t > zero, "t > 0 (infinitesimal positive)");
    check(!(t < zero), "not t < 0");
    check((-t) < zero, "-t < 0");
    J<T, Nd> f = t - t * t;
    check(f > zero, "t - t^2 > 0");
    J<T, Nd> g = -t + J<T, Nd>(5) * t * t;
    check(g < zero, "-t + 5 t^2 < 0 (leading term wins)");
    check(g.sign() == -1, "sign(-t + 5 t^2) == -1");
    check(J<T, Nd>(3) > J<T, Nd>(2), "constant 3 > 2");
    check(f == (t - t * t), "equality of identical jets");
    check(zero >= zero, "0 >= 0");
    check(f >= zero, "t - t^2 >= 0");
  }

  // Matrix determinant / derivatives over jets.
  {
    constexpr int Nd = 3;
    J<T, Nd> t = J<T, Nd>::var();
    J<T, Nd> m00 = J<T, Nd>(1) + t, m01 = t, m10 = t, m11 = J<T, Nd>(1) - t;
    J<T, Nd> det = m00 * m11 - m01 * m10; // 1 - 2 t^2
    check(det[0] == T(1) && det[1] == T(0) && det[2] == T(-2) && det[3] == T(0),
          "2x2 jet determinant 1 - 2 t^2");
    check(jet_deriv(det, 0) == T(1) && jet_deriv(det, 1) == T(0) &&
              jet_deriv(det, 2) == T(-4),
          "jet_deriv of 1 - 2 t^2 (0th, 1st, 2nd)");
  }

  // MyMatrix<jet> with DeterminantMat / Inverse (Q + t H, Q = [[2,-1],[-1,2]]).
  {
    constexpr int Nd = 2;
    using JT = jet<T, Nd>;
    JT t = JT::var();
    MyMatrix<JT> M(2, 2);
    M(0, 0) = JT(2);
    M(0, 1) = JT(-1) + t;
    M(1, 0) = JT(-1) + t;
    M(1, 1) = JT(2);
    JT det = DeterminantMat(M); // 3 + 2t - t^2
    check(det[0] == T(3) && det[1] == T(2) && det[2] == T(-1),
          "DeterminantMat over jets: 3 + 2t - t^2");
    MyMatrix<JT> Inv = Inverse(M);
    MyMatrix<JT> Prod = Inv * M;
    check(Prod(0, 0) == JT(1) && Prod(1, 1) == JT(1) && Prod(0, 1) == JT(0) &&
              Prod(1, 0) == JT(0),
          "Inverse over jets: Inv * M == I");
  }

  // The pivot-selection fix. M = [[t, 1], [1, t]]: the (0,0) entry is t (order 1,
  // NOT invertible), the (0,1) entry is 1 (order 0). The historical
  // "first non-zero pivot" rule would pick column 0 = t and try 1/t (a Laurent
  // series -> failure); the jet "smallest order" rule picks column 1. det = t^2
  // - 1 = -1 + t^2, invertible, so this must succeed and be exact.
  {
    constexpr int Nd = 2;
    using JT = jet<T, Nd>;
    JT t = JT::var();
    MyMatrix<JT> M(2, 2);
    M(0, 0) = t;
    M(0, 1) = JT(1);
    M(1, 0) = JT(1);
    M(1, 1) = t;
    JT det = DeterminantMat(M); // -1 + t^2
    check(det[0] == T(-1) && det[1] == T(0) && det[2] == T(1),
          "pivot fix: det [[t,1],[1,t]] = -1 + t^2");
    MyMatrix<JT> Inv = Inverse(M);
    MyMatrix<JT> Prod = Inv * M;
    check(Prod(0, 0) == JT(1) && Prod(1, 1) == JT(1) && Prod(0, 1) == JT(0) &&
              Prod(1, 0) == JT(0),
          "pivot fix: Inverse [[t,1],[1,t]] * M == I");
    // SolutionMat: solve x M = (1, 0) (exact rational function of t).
    MyVector<JT> rhs(2);
    rhs(0) = JT(1);
    rhs(1) = JT(0);
    std::optional<MyVector<JT>> sol = SolutionMat(M, rhs);
    bool sol_ok = sol.has_value();
    if (sol_ok) {
      // x = rhs * M^{-1}; verify x M = rhs.
      MyVector<JT> xM = ProductVectorMatrix(*sol, M);
      sol_ok = (xM(0) == JT(1)) && (xM(1) == JT(0));
    }
    check(sol_ok, "pivot fix: SolutionMat over jets on [[t,1],[1,t]]");
  }

  // Division-free fallback: a jet matrix whose DETERMINANT itself has c0 == 0
  // (its constant-term matrix is singular). Elimination would be forced onto a
  // zero-divisor pivot, so DeterminantMat must detect this and switch to the
  // division-free permutation expansion. This is the "split simplex" case
  // det = alpha t^d + ... that vanishes at t = 0 but not for t > 0.
  {
    constexpr int Nd = 2;
    using JT = jet<T, Nd>;
    JT t = JT::var();
    // [[t, 0], [0, t]] -> det = t^2 (valuation 2; constant-term matrix all zero).
    MyMatrix<JT> M(2, 2);
    M(0, 0) = t;
    M(0, 1) = JT(0);
    M(1, 0) = JT(0);
    M(1, 1) = t;
    JT det = DeterminantMat(M);
    check(det[0] == T(0) && det[1] == T(0) && det[2] == T(1),
          "division-free det [[t,0],[0,t]] = t^2");
    // [[t, 1], [t, t]] -> det = t^2 - t (valuation 1).
    MyMatrix<JT> M2(2, 2);
    M2(0, 0) = t;
    M2(0, 1) = JT(1);
    M2(1, 0) = t;
    M2(1, 1) = t;
    JT det2 = DeterminantMat(M2);
    check(det2[0] == T(0) && det2[1] == T(-1) && det2[2] == T(1),
          "division-free det [[t,1],[t,t]] = -t + t^2");
  }

  // Berkowitz cross-check: on ordinary rationals it must equal the elimination
  // determinant, for several sizes.
  {
    bool ok = true;
    auto fill = [](int nn, int seed) {
      MyMatrix<T> A(nn, nn);
      for (int i = 0; i < nn; i++)
        for (int j = 0; j < nn; j++)
          A(i, j) = T(((i * 7 + j * 3 + seed * 5) % 11) - 5);
      return A;
    };
    for (int nn = 1; nn <= 6 && ok; nn++)
      for (int seed = 0; seed < 4 && ok; seed++) {
        MyMatrix<T> A = fill(nn, seed);
        if (DeterminantMatBerkowitz(A) != DeterminantMatKernel(A))
          ok = false;
      }
    check(ok, "Berkowitz == elimination determinant (n=1..6)");
  }

  // Higher-corank jet case: a 3x3 that is corank 3 at t = 0 (whole matrix is
  // O(t)), exercising the full Berkowitz residual. M = t*I3 -> det = t^3.
  {
    constexpr int Nd = 4;
    using JT = jet<T, Nd>;
    JT t = JT::var();
    MyMatrix<JT> M(3, 3);
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++)
        M(i, j) = (i == j) ? t : JT(0);
    JT det = DeterminantMat(M);
    check(det[0] == T(0) && det[1] == T(0) && det[2] == T(0) && det[3] == T(1),
          "division-free det t*I3 = t^3 (corank-3 Berkowitz residual)");
  }

  if (n_error == 0)
    std::cerr << "All jet_number tests passed.\n";
  else
    std::cerr << n_error << " jet_number test(s) FAILED.\n";
  return n_error == 0 ? 0 : 1;
}
