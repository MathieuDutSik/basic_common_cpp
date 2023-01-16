// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Heuristic_ThompsonSampling.h"

int main() {
  size_t n_max = 10;
  LimitedEmpiricalDistributionFunction ledf(n_max);
  size_t N_split = 4587;
  for (int i = 0; i < 100; i++) {
    double new_val = static_cast<double>(10) * static_cast<double>(random() % N_split) / static_cast<double>(N_split);
    std::cerr << "i=" << i << " new_val=" << new_val << "\n";
    ledf.insert_value(new_val);
  }
  size_t N_test = 17;
  for (int u = 0; u < 100; u++) {
    double alpha = static_cast<double>(random() % N_test) / static_cast<double>(N_test);
    double samp = ledf.get_percentile(alpha);
    std::cerr << "u=" << u << " alpha=" << alpha << " samp=" << samp << "\n";
  }
  std::cerr << "ledf=" << ledf.string() << "\n";
}
