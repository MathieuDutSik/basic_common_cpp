// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "NumberGivaro.h"

int main() {
  for (int i = 0; i < 100; i++) {
    int eval_i = random() % 3 - 1;
    mpz_class eval = eval_i;
    for (int k = 0; k < 10; k++)
      eval *= random();
    Givaro::Integer eval_g = GetGivaroInteger(eval);
    mpz_class evalb = ConvertGivaroInteger(eval_g);
    std::cerr << "eval=" << eval << "\n";
    if (eval != evalb) {
      std::cerr << "i=" << i << " eval=" << eval << " evalb=" << evalb << "\n";
      std::cerr << "Error in the system\n";
      return 44;
    }
  }
  std::cerr << "Normal termination of the program\n";
}
