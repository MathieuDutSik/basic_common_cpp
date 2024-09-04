// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalFunctions.h"
#include "GRAPH_traces.h"
#include "Temp_common.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "GRAPH_Kconnectivity [FileGraph] [k]\n";
      std::cerr << "\n";
      return -1;
    }
    using Tgr = GraphBitset;
    //
    std::string FileGraph = argv[1];
    Tgr eGR = GRAPH_ReadFile<GraphBitset>(FileGraph);
    GRAPH_Write(std::cerr, eGR);
    //
    std::string k_str = argv[2];
    size_t k = ParseScalar<size_t>(k_str);
    //
    bool result = IsKConnectedGraph(eGR, k);
    std::cerr << "result=" << result << "\n";
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something went wrong in the code\n";
    exit(e.eVal);
  }
}
