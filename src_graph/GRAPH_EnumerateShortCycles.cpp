#include "Basic_file.h"
#include "GRAPH_GraphicalFunctions.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "This program is used as\n";
      std::cerr << "GRAPH_EnumerateShortCycles [DataGraph]\n";
      std::cerr << "\n";
      std::cerr << "DataGraph : The file containing the graph\n";
      return -1;
    }
    //
    std::ifstream GRAfs(argv[1]);
    GraphBitset eGR=GRAPH_Read<GraphBitset>(GRAfs);
    GRAPH_PrintOutput(std::cerr, eGR);
    //
    std::vector<std::vector<size_t>> ListCycles = GRAPH_FindAllCycles(eGR);
    size_t nbCycle=ListCycles.size();
    std::cerr << "nbCycle=" << nbCycle << "\n";
    for (size_t i_cyc=0; i_cyc<nbCycle; i_cyc++) {
      std::vector<size_t> const& eCyc = ListCycles[i_cyc];
      std::cerr << "i_cyc=" << i_cyc << " eCyc =";
      for (auto & eVal : eCyc)
        std::cerr << " " << eVal;
      std::cerr << "\n";
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
