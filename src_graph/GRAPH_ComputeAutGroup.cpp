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
      std::cerr << "GRAPH_ComputeAutGroup [DataGraph] [opt]\n";
      std::cerr << "\n";
      return -1;
    }
    std::cerr << "Reading input\n";
    //
    std::ifstream GRAfs(argv[1]);
    int opt;
    sscanf(argv[2], "%d", &opt);
    //
    GraphBitset eGR=GRAPH_Read<GraphBitset>(GRAfs);
    std::vector<std::vector<unsigned int>> ListGen;
    if (opt == 1) {
      std::cerr << "Running TRACES GetCanonicalOrdering\n";
      ListGen = TRACES_GetListGenerators(eGR);
    } else {
      std::cerr << "Running BLISS GetCanonicalOrdering\n";
      ListGen = BLISS_GetListGenerators(eGR);
    }
    int nbVert=eGR.GetNbVert();
    std::cout << "local ListGen;\n";
    std::cout << "ListGen:=[";
    bool IsFirst=true;
    for (auto & eGen : ListGen) {
      if (!IsFirst)
        std::cout << ",\n";
      IsFirst=false;
      std::cout << "[";
      for (int i=0; i<nbVert; i++) {
        if (i>0)
          std::cout << ",";
        int eVal = 1 + eGen[i];
        std::cout << eVal;
      }
      std::cout << "]";
    }
    std::cout << "];\n";
    std::cout << "return Group(List(ListGen,PermList));\n";
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
