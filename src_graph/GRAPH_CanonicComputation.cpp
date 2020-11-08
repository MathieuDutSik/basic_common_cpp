#include "GRAPH_nauty.h"

int main(int argc, char *argv[])
{
  try {
    if (argc != 2) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "GRAPH_CanonicOrdering [DataGraph]\n";
      std::cerr << "\n";
      return -1;
    }
    std::cerr << "Reading input\n";
    //
    std::ifstream GRAfs(argv[1]);
    GraphBitset eGR=GRAPH_Read<GraphBitset>(GRAfs);
    std::vector<int> V = GetNautyCanonicalOrdering(eGR);
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
