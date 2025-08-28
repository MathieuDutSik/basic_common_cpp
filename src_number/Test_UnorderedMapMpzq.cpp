// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "Temp_common.h"
// clang-format on

int main() {
  int n = 1000;
  int expo = 10;
  std::unordered_map<mpz_class, int> M1;
  std::unordered_map<mpq_class, int> M2;
  auto get_mpz = [&]() -> mpz_class {
    int eval_i = random() % 3 - 1;
    mpz_class eval = eval_i;
    for (int k = 0; k < expo; k++)
      eval *= (1 + random() % 100);
    return eval;
  };
  auto get_mpq = [&]() -> mpq_class {
    mpq_class val1 = mpq_class(get_mpz());
    mpq_class val2 = mpq_class(get_mpz());
    if (val2 > 0)
      return val1 / val2;
    return val1 / (1 - val2);
  };
  std::set<mpz_class> S1;
  std::set<mpq_class> S2;
  for (int i = 0; i < n; i++) {
    mpz_class val1 = get_mpz();
    int int1 = random() / 100;
    M1[val1] = int1;
    S1.insert(val1);
    //
    mpq_class val2 = get_mpq();
    int int2 = random() / 100;
    M2[val2] = int2;
    S2.insert(val2);
  }
  std::cerr << "|M1|=" << M1.size() << " |M2|=" << M2.size()
            << " |S1|=" << S1.size() << " |S2|=" << S2.size() << "\n";
  std::cerr << "Normal termination of the program\n";
}
