// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>

#include "Heuristic_ThompsonSampling.h"
#include "NumberTheoryGmp.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time;
  try {
    //  using T = mpq_class;
    using T = mpz_class;
    //    using T = T_uint64_t;
    // The chosenoptions
    FullNamelist eFull = NAMELIST_ThompsonSamplingRuntime();
    if (argc != 2) {
      std::cerr << "Program is used as\n";
      std::cerr << "Test_Thompson_Sampling [FileI]\n";
      std::cerr << "File format is\n";
      eFull.NAMELIST_WriteNamelistFile(std::cerr, true);
      throw TerminalException{1};
    }
    std::string filename = argv[1];
    NAMELIST_ReadNamelistFile(filename, eFull);

    ThompsonSamplingHeuristic<T> TSH(eFull, std::cerr);

    std::vector<std::string> l_input = GetHeuristicInput(TSH.heu);
    std::map<std::string, std::vector<T>> l_poss;
    for (auto &eKey : l_input)
      l_poss[eKey] = GetHeuristicPivots(TSH.heu, eKey);
    size_t N = 10000;
    std::map<std::string, T> TheCand;
    for (size_t i = 0; i < N; i++) {
      for (auto &eKey : l_input) {
        auto &e_vect = l_poss[eKey];
        size_t len = e_vect.size();
        size_t pos = random() % len;
        T val = e_vect[pos];
        TheCand[eKey] = val;
      }
      std::string choice = TSH.get_eval(TheCand);
      TSH.pop(std::cerr);
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
