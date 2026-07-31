// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Fp.h"

const long P = 2147389441;

int main(int argc, char *argv[]) {
  using Tint = Fp<long, P>;

  Tint x = Tint(5);
  Tint y = Tint(-3 * P + 5);
  std::cerr << (x == y) << std::endl;

  std::cerr << (x + y - x == y) << std::endl;

  Tint xinv = Tint(1) / x;
  std::cerr << (x.get_num() * xinv.get_num()) % P << std::endl;

  std::cerr << (Tint(7) * xinv == Tint(7) / x) << std::endl;

  std::pair<long, long> r = (Tint(7) / Tint(5)).rational_lift();
  std::cout << (r.first * 5 == r.second * 7) << std::endl;

  std::pair<long, long> r2 = Tint(2131).rational_lift();
  std::cerr << r2.first << " / " << r2.second << std::endl;

  return 0;
}
