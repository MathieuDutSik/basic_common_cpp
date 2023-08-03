// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Boost_bitset_kernel.h"
#include "MAT_Matrix.h"

int main() {
  size_t n = 20;
  Face f(n);
  for (size_t i = 0; i < n; i++)
    f[i] = random() % 2;

  std::cerr << "f (V1) =";
  for (size_t i = 0; i < n; i++)
    std::cerr << " " << f[i];
  std::cerr << "\n";

  std::cerr << "f (V2) =";
  for (size_t i = 0; i < n; i++)
    if (f[i] == 1)
      std::cerr << " " << i;
  std::cerr << "\n";

  //  std::cerr << "f (V3) =";
  //  for (auto& eVal : f) {
  //    std::cerr << " " << eVal;
  //  }
  //  std::cerr << "\n";
}
