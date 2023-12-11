// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Heuristic_ThompsonSampling.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time1;
  try {
    if (argc != 2) {
      std::cerr << "Analysis_EmpiricalDistribution [FileI]\n";
      throw TerminalException{1};
    }
    std::string FileI = argv[1];
    std::vector<std::string> ListLines = ReadFullFile(FileI);
    for (auto &eLine : ListLines) {
      std::vector<std::string> LStr = STRING_Split(eLine, "ListValWei = ");
      if (LStr.size() == 2) {
        std::string desc = LStr[1];
        LimitedEmpiricalDistributionFunction ledf(25, 0, "sampled", desc);
        double average = ledf.get_average();
        std::cerr << "average=" << average << "\n";
      }
    }
    std::cerr << "Normal termination of Analysis_EmpiricalDistribution\n";
  } catch (TerminalException const &e) {
    std::cerr << "Error in Analysis_EmpiricalDistribution\n";
    exit(e.eVal);
  }
  runtime(time1);
}
