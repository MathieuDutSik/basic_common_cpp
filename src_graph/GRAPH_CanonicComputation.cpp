#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalFunctions.h"
#include "GRAPH_bliss.h"
#include "GRAPH_traces.h"
#include "Temp_common.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "GRAPH_CanonicOrdering [DataGraph] [opt]\n";
      std::cerr << "\n";
      return -1;
    }
    using Tgr = GraphBitset;
    using Tidx = unsigned int;
    //
    std::ifstream GRAfs(argv[1]);
    Tgr eGR = GRAPH_Read<Tgr>(GRAfs);
    GRAPH_Write(std::cerr, eGR);
    //
    int opt;
    sscanf(argv[2], "%d", &opt);
    std::cerr << "opt=" << opt << "\n";
    //
    std::vector<Tidx> V;
    if (opt == 1) {
      std::cerr << "Running TRACES GetCanonicalOrdering\n";
      V = TRACES_GetCanonicalOrdering<Tgr, Tidx>(eGR);
    }
    if (opt == 2) {
      std::cerr << "Running TRACES GetCanonicalOrdering_Arr_Test\n";
      V = TRACES_GetCanonicalOrdering_Arr_Test<Tgr, Tidx>(eGR);
    }
    if (opt == 3) {
      std::cerr << "Running BLISS GetCanonicalOrdering\n";
      V = BLISS_GetCanonicalOrdering<Tgr, Tidx>(eGR);
    }
    size_t nbVert = V.size();
    std::cout << "return [";
    for (size_t i = 0; i < nbVert; i++) {
      if (i > 0)
        std::cout << ",";
      std::cout << V[i];
    }
    std::cout << "];\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
