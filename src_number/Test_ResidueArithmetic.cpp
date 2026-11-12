// Copyright (C) 2026 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "ResidueArithmetic.h"
#include <random>
// clang-format on

/*
  The residue operations of ResidueArithmetic.h against exact integers, on
  random residues for moduli up to 2^62, where every naive product would
  overflow, and the refusal of a modulus beyond 2^62.
 */

int main() {
  try {
    using Tint = mpz_class;
    std::mt19937_64 rng(1);
    std::vector<int64_t> moduli{1, 2, 3, 4, 9, 1024, 3486784401,
                                (int64_t(1) << 31) + 11, int64_t(3486784401) * 3,
                                kMaxResidueModulus - 1, kMaxResidueModulus,
                                ResiduePrimePower(3, 39), ResiduePrimePower(2, 62)};
    size_t n_check = 0;
    for (int64_t m : moduli) {
      Tint m_T = UniversalScalarConversion<Tint, int64_t>(m);
      std::uniform_int_distribution<int64_t> unif(0, m - 1);
      for (int iter = 0; iter < 2000; iter++) {
        int64_t a = unif(rng), b = unif(rng);
        Tint a_T = UniversalScalarConversion<Tint, int64_t>(a);
        Tint b_T = UniversalScalarConversion<Tint, int64_t>(b);
        auto exact = [&](Tint const &x) -> int64_t { return ResModOf(x, m); };
        if (ResMul(a, b, m) != exact(a_T * b_T) ||
            ResAdd(a, b, m) != exact(a_T + b_T) ||
            ResSub(a, b, m) != exact(a_T - b_T) ||
            ResNeg(a, m) != exact(-a_T)) {
          std::cerr << "Test_ResidueArithmetic: wrong result for a=" << a
                    << " b=" << b << " m=" << m << "\n";
          return 1;
        }
        Tint g;
        mpz_gcd(g.get_mpz_t(), a_T.get_mpz_t(), m_T.get_mpz_t());
        if (m > 1 && g == 1) {
          int64_t ai = ResInv(a, m);
          if (exact(a_T * UniversalScalarConversion<Tint, int64_t>(ai)) != 1) {
            std::cerr << "Test_ResidueArithmetic: wrong inverse of " << a
                      << " modulo " << m << "\n";
            return 1;
          }
        }
        n_check++;
      }
    }
    // a modulus beyond 2^62 is refused
    bool refused = false;
    try {
      (void)ResiduePrimePower(2, 63);
    } catch (TerminalException const &) {
      refused = true;
    }
    try {
      Tint big = UniversalScalarConversion<Tint, int64_t>(kMaxResidueModulus);
      big += 1;
      (void)ResidueModulusOf(big);
      refused = false;
    } catch (TerminalException const &) {
    }
    if (!refused) {
      std::cerr << "Test_ResidueArithmetic: a modulus beyond 2^62 was accepted\n";
      return 1;
    }
    std::cerr << "Test_ResidueArithmetic: " << n_check
              << " random checks passed, moduli beyond 2^62 refused\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Test_ResidueArithmetic\n";
    exit(e.eVal);
  }
}
