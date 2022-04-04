#include "Coloring.h"
#include "GRAPH_BitsetType.h"
#include "GRAPH_GraphicalFunctions.h"
#include "MAT_Matrix.h"

int main(int argc, char *argv[]) {
  try {
    if (argc > 4 || argc < 2) {
      std::cerr << "compute_chromatic_number is used as\n";
      std::cerr
          << "compute_chromatic_number [GRA] [StartChromatic] [OutFile]\n";
      std::cerr << "or\n";
      std::cerr << "compute_chromatic_number [GRA] [StartChromatic]\n";
      std::cerr << "or\n";
      std::cerr << "compute_chromatic_number [GRA]\n";
      throw TerminalException{1};
    }
    //
    std::string FileName = argv[1];
    int StartChromatic = 1;
    if (argc >= 3) {
      sscanf(argv[2], "%d", &StartChromatic);
    }
    std::cerr << "Starting chromatic number is " << StartChromatic << "\n";
    //
    std::ifstream GRAfs(FileName);
    GraphBitset eGR = GRAPH_Read<GraphBitset>(GRAfs);

    int nbColor = StartChromatic;
    while (true) {
      std::pair<bool, std::vector<int>> ePair =
          GetColoringOrFail<GraphBitset>(eGR, nbColor);
      if (ePair.first) {
        std::vector<int> V = ePair.second;
        int nbPoint = V.size();
        std::cerr << "The Chromatic number is " << nbColor << "\n";
        if (argc == 4) {
          std::string OutFile = argv[3];
          std::ofstream os(OutFile);
          os << nbColor << "\n";
          os << nbPoint << "\n";
          for (int iPoint = 0; iPoint < nbPoint; iPoint++)
            os << V[iPoint] << "\n";
        }
        break;
      }
      nbColor++;
    }

  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
