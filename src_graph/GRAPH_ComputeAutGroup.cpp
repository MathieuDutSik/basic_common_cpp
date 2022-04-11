#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalFunctions.h"
#include "GRAPH_bliss.h"
#include "GRAPH_traces.h"
#include "Temp_common.h"

int main(int argc, char *argv[]) {
  try {
    if (argc != 3) {
      std::cerr << "Number of argument is = " << argc << "\n";
      std::cerr << "GRAPH_ComputeAutGroup [DataGraph] [opt]\n";
      std::cerr << "\n";
      return -1;
    }
    using Tgr = GraphBitset;
    using Tidx = unsigned int;
    //
    std::ifstream GRAfs(argv[1]);
    Tgr eGR = GRAPH_Read<GraphBitset>(GRAfs);
    GRAPH_Write(std::cerr, eGR);
    //
    int opt;
    sscanf(argv[2], "%d", &opt);
    std::cerr << "opt=" << opt << "\n";
    //
    size_t nbRow = eGR.GetNbVert();
    std::vector<std::vector<Tidx>> ListGen;
    if (opt == 1) {
      std::cerr << "Running TRACES_GetListGenerators\n";
      ListGen = TRACES_GetListGenerators<Tgr, Tidx>(eGR, nbRow);
    }
    if (opt == 2) {
      std::cerr << "Running TRACES_GetListGenerators\n";
      ListGen = TRACES_GetListGenerators_Arr_Test<Tgr, Tidx>(eGR, nbRow);
    }
    if (opt == 3) {
      std::cerr << "Running BLISS_GetListGenerators\n";
      ListGen = BLISS_GetListGenerators<Tgr, Tidx>(eGR, nbRow);
    }
    size_t nbVert = eGR.GetNbVert();
    std::cout << "local ListGen;\n";
    std::cout << "ListGen:=[";
    bool IsFirst = true;
    for (auto &eGen : ListGen) {
      if (!IsFirst)
        std::cout << ",\n";
      IsFirst = false;
      std::cout << "[";
      for (size_t i = 0; i < nbVert; i++) {
        if (i > 0)
          std::cout << ",";
        int eVal = 1 + int(eGen[i]);
        std::cout << eVal;
      }
      std::cout << "]";
    }
    std::cout << "];\n";
    std::cout << "return Group(List(ListGen,PermList));\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
