#include "NumberTheory.h"
#include "factorizations.h"

int main(int argc, char *argv[]) {
  using T = mpz_class;
  // using T=mpq_class;
  try {
    for (int n = 1; n < 500; n++) {
      T n_T = n;
      std::vector<T> V = FactorsInt(n_T);
      std::cerr << "n=" << n_T << " Fact=";
      for (auto &val : V)
        std::cerr << " " << val;
      std::cerr << "\n";
      std::vector<T> Ldiv = GetAllFactors(n_T);
      std::cerr << "  divisors =";
      for (auto &val : Ldiv)
        std::cerr << " " << val;
      std::cerr << "\n";
    }
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
