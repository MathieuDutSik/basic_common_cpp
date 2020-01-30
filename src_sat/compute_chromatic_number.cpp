#include "Coloring.h"
#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalFunctions.h"
#include "MAT_Matrix.h"

int main(int argc, char* argv[])
{
  try {
    if (argc != 4) {
      std::cerr << "compute_chromatic_number is used as\n";
      std::cerr << "compute_chromatic_number [GRA] [StartChromatic] [OutFile]\n";
      throw TerminalException{1};
    }
    //
    std::string FileName = argv[1];
    int StartChromatic;
    sscanf(argv[2], "%d", &StartChromatic);
    std::string OutFile = argv[3];
    //
    std::ifstream GRAfs(FileName);
    GraphBitset eGR=GRAPH_Read<GraphBitset>(GRAfs);

    int nbColor = StartChromatic;
    while(true) {
      std::pair<bool, std::vector<int>> ePair = GetColoringOrFail<GraphBitset>(eGR, nbColor);
      if (ePair.first) {
        std::vector<int> V = ePair.second;
        int nbPoint=V.size();
        std::ofstream os(OutFile);
        os << nbColor << "\n";
        os << nbPoint << "\n";
        for (int iPoint=0; iPoint<nbPoint; iPoint++)
          os << V[iPoint] << "\n";
        break;
      }
      nbColor++;
    }

  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
