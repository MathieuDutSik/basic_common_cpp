// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryPadic.h"
#include "TypeConversion.h"
// clang-format on

template <typename T> void process() {
  std::vector<T> ListVal;
  int n_vert=100;
  for (int i=0; i<n_vert; i++) {
    int val = rand() % 100000;
    T val_T = UniversalScalarConversion<T,int>(val);
    ListVal.push_back(val_T);
  }
  std::vector<int> primes{2, 3, 5, 7, 23, 719, 6991};
  for (auto & p_i : primes) {
    T p(p_i);
    std::vector<Padic<T>> ListVal_adic;
    for (auto & eVal1 : ListVal) {
      Padic<T> eVal2 = Padic_from_integer(eVal1, p);
      Padic<T> eVal3 = Padic_reduce_precision(eVal2, 3);
      ListVal_adic.push_back(eVal3);
    }
    for (int iVert=0; iVert<n_vert; iVert++) {
      for (int jVert=iVert+1; jVert<n_vert; jVert++) {
        Padic<T> ValJ_inv = Padic_inverse(ListVal_adic[jVert], p);
        Padic<T> prod = Padic_product(
      }
    }
  }
}

int main(int argc, char *argv[]) {
  try {
    process<mpz_class>();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "There is a bug to be resolved in the code\n";
    exit(e.eVal);
  }
}
