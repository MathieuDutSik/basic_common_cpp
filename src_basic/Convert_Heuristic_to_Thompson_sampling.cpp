// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "Heuristic_fct.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "Convert_Heuristic_to_Thompson_sampling [fileI] [FileO]\n";
      throw TerminalException{1};
    }
    using T = T_uint64_t;
    std::string FileI = argv[1];
    std::string FileO = argv[2];
    //
    std::ifstream is(FileI);
    TheHeuristic<T> heu = ReadHeuristic<T>(is);
    FullNamelist eFull = ConvertHeuristicToFullNamelist(heu);
    //
    std::ofstream os(FileO);
    NAMELIST_WriteNamelistFile(os, eFull);
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Erroneous termination of the program\n";
    exit(e.eVal);
  }
}
