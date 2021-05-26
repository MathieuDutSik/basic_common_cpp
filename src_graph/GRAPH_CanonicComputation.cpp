#include "GRAPH_traces.h"
#include "GRAPH_bliss.h"
#include "Temp_common.h"
#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalFunctions.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "GRAPH_CanonicOrdering [DataGraph] [opt]\n";
      std::cerr << "\n";
      return -1;
    }
    std::cerr << "Reading input\n";
    //
    std::ifstream GRAfs(argv[1]);
    int opt;
    sscanf(argv[2], "%d", &opt);
    //
    using Tgr = GraphBitset;
    using Tidx = unsigned int;
    Tgr eGR=GRAPH_Read<Tgr>(GRAfs);
    std::vector<Tidx> V;
    if (opt == 1) {
      std::cerr << "Running TRACES GetCanonicalOrdering\n";
      V = TRACES_GetCanonicalOrdering<Tgr,Tidx>(eGR);
    } else {
      std::cerr << "Running BLISS GetCanonicalOrdering\n";
      V = BLISS_GetCanonicalOrdering<Tgr,Tidx>(eGR);
    }
    int nbVert=V.size();
    std::cout << "return [";
    for (int i=0; i<nbVert; i++) {
      if (i>0)
        std::cout << ",";
      std::cout << V[i];
    }
    std::cout << "];\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
