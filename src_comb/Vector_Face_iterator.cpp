// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Boost_bitset.h"

int main() {
  size_t n = 10;
  size_t n_iter = 20;
  vectface vf(n);

  for (size_t iter = 0; iter < n_iter; iter++) {
    Face f(n);
    for (size_t i = 0; i < n; i++) {
      f[i] = random() % 2;
    }
    vf.push_back(f);
  }
  //
  //  std::sort(vf.begin(), vf.end());

  for (auto &f : vf) {
    std::cerr << " |f|=" << f.size() << " / " << f.count() << "\n";
  }
}
