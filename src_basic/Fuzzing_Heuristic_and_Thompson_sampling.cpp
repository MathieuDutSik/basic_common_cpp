// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Heuristic_ThompsonSampling.h"

int main(int argc, char *argv[]) {
  srand_random_set();
  HumanTime time;
  try {
    if (argc != 2) {
      std::cerr << "Fuzzing_Heuristic_and_Thompson_sampling [fileI]\n";
      throw TerminalException{1};
    }
    using T = T_uint64_t;
    std::string FileI = argv[1];
    //
    std::ifstream is(FileI);
    TheHeuristic<T> heu = ReadHeuristic<T>(is);
    FullNamelist eFull = ConvertHeuristicToFullNamelist(heu);
    std::cerr << "heu=\n" << heu << "\n";
    std::cerr << "eFull=\n";
    NAMELIST_WriteNamelistFile(std::cerr, eFull, false);
    //
    ThompsonSamplingHeuristic<T> TSH(eFull, std::cerr);
    std::cerr << "The TSH object has been BUILD\n";
    std::vector<std::string> l_input = GetHeuristicInput(TSH.heu);
    std::map<std::string, std::vector<T>> l_poss;
    for (auto &eKey : l_input)
      l_poss[eKey] = GetHeuristicPivots(TSH.heu, eKey);
    //
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
      std::string choice_TS = TSH.get_eval(TheCand);
      TSH.pop(std::cerr);
      std::string choice_Heu = HeuristicEvaluation(TheCand, heu);
      if (choice_TS != choice_Heu) {
        std::cerr << "The heuristic returned different values from the "
                     "heuristic and TS. Bug to resolve\n";
        throw TerminalException{1};
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
  runtime(time);
}
