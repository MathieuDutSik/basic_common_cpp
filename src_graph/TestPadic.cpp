// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryPadic.h"
#include "TypeConversion.h"
#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalFunctions.h"
// clang-format on

template <typename T> void process() {
  std::vector<T> ListVal;
  int n_vert = 100;
  for (int i = 0; i < n_vert; i++) {
    int val = rand() % 1000;
    std::cerr << "i=" << i << " val=" << val << "\n";
    T val_T = UniversalScalarConversion<T, int>(val);
    ListVal.push_back(val_T);
  }
  std::vector<int> primes{2, 3, 5, 7, 23, 719, 6991};
  // std::vector<int> primes{2};
  for (auto &p_i : primes) {
    T p(p_i);
    std::cerr << "-------------------\n";
    std::cerr << "p=" << p << "\n";
    std::vector<Padic<T>> ListVal_adic;
    for (auto &eVal1 : ListVal) {
      std::cerr << "eVal1=" << eVal1 << "\n";
      Padic<T> eVal2 = Padic_from_positive_integer(eVal1, p);
      std::cerr << "  eVal2 : ";
      Padic_debug_print(eVal2, std::cerr);
      Padic<T> eVal3 = Padic_reduce_precision(eVal2, 3);
      std::cerr << "  eVal3 : ";
      Padic_debug_print(eVal3, std::cerr);
      ListVal_adic.push_back(eVal3);
    }
    using Tgr = GraphBitset;
    Tgr eGR(n_vert);
    for (int iVert = 0; iVert < n_vert; iVert++) {
      for (int jVert = iVert + 1; jVert < n_vert; jVert++) {
        std::cerr << "iVert=" << iVert << " jVert=" << jVert << "\n";
        std::cerr << "  V[jVert] : ";
        Padic_debug_print(ListVal_adic[jVert], std::cerr);
        std::cerr << "  V[iVert] : ";
        Padic_debug_print(ListVal_adic[iVert], std::cerr);
        Padic<T> ValJ_inv = Padic_inverse(ListVal_adic[jVert], p);
        std::cerr << "  ValJ_inv : ";
        Padic_debug_print(ValJ_inv, std::cerr);
        Padic<T> prod = Padic_product(ListVal_adic[iVert], ValJ_inv, p);
        std::cerr << "  prod : ";
        Padic_debug_print(prod, std::cerr);
        bool test = Padic_is_square(prod, p);
        std::cerr << "  test=" << test << "\n";
        if (test) {
          eGR.AddAdjacent(iVert, jVert);
          eGR.AddAdjacent(jVert, iVert);
        }
      }
    }
    std::vector<std::vector<size_t>> LConn = ConnectedComponents_set(eGR);
    std::cerr << "p=" << p << " |LConn|=" << LConn.size() << "\n";
    for (auto &eConn : LConn) {
      std::cerr << " " << eConn.size();
    }
    std::cerr << "\n";
  }
}

int main() {
  try {
    process<mpz_class>();
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "There is a bug to be resolved in the code\n";
    exit(e.eVal);
  }
}
