// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryPadic.h"
#include "TypeConversion.h"
// clang-format on

int main() {
  using T = mpz_class;
  try {
    for (int p_i=2; p_i<30; p_i++) {
      T p(p_i);
      std::vector<T> classes = Padic_get_residue_classes(p);
      auto get_class=[&](Padic<T> const& u) -> size_t {
        Padic<T> u_red = Padic_reduce_precision(u, 3);
        for (auto & cand : classes) {
          Padic<T> v = Padic_from_positive_integer(cand, p);
          Padic<T> v_red = Padic_reduce_precision(v, 3);
          Padic<T> prod = Padic_product(v_red, u_red, p);
          bool test = Padic_is_square(prod, p);
          if (test) {
            std::pair<size_t, T> pair = Padic_decompose(cand, p);
            return pair.first;
          }
        }
        std::cerr << "Failed to find a class\n";
        throw TerminalException{1};
      };
      if (IsPrime(p)) {
        int val_i = 1 + rand() % 1000;
        T val1(val_i);
        Padic<T> val2 = Padic_from_positive_integer(val1, p);
        size_t pos = get_class(val2);
        std::cerr << "pos=" << pos << "\n";
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the computation\n";
    exit(e.eVal);
  }
}
