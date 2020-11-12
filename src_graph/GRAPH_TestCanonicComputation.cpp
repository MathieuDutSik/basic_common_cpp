#include "GRAPH_nauty.h"
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
    GraphBitset eGR=GRAPH_Read<GraphBitset>(GRAfs);
    //
    auto get_string_expression=[&](GraphBitset const& eGR) -> std::string {
      std::vector<int> V;
      if (opt == 1) {
        V = NAUTY_GetCanonicalOrdering(eGR);
      } else {
        V = BLISS_GetCanonicalOrdering(eGR);
      }
      int nof_vertices = V.size();
      std::vector<int> clR(nof_vertices);
      for (int i=0; i<nof_vertices; i++)
        clR[V[i]]=i;
      std::string strRet;
      for (int iVert=0; iVert<nof_vertices; iVert++) {
        int iVertCan = clR[iVert];
        strRet += " (";
        if (eGR.GetHasVertexColor()) {
          int eColor = eGR.GetColor(iVertCan);
          strRet += std::to_string(eColor) + ":";
        }
        //
        for (int jVert=0; jVert<nof_vertices; jVert++) {
          int jVertCan = clR[jVert];
          bool eVal_b = eGR.IsAdjacent(iVertCan, jVertCan);
          strRet += " " + std::to_string(eVal_b);
        }
        strRet += ")";
      }
      return strRet;
    };
    auto random_perm_graph=[&](GraphBitset const& eGR) -> GraphBitset {
      int nbVert = eGR.GetNbVert();
      std::vector<int> ePerm = RandomPermutation(nbVert);
      GraphBitset eGR2(nbVert);
      for (int iVert=0; iVert<nbVert; iVert++) {
        int iVert2 = ePerm[iVert];
        for (auto& eAdj : eGR.Adjacency(iVert)) {
          int eAdj2 = ePerm[eAdj];
          eGR2.AddAdjacent(iVert2, eAdj2);
        }
      }
      return eGR2;
    };
    std::string str1 = get_string_expression(eGR);
    for (int i=0; i<10; i++) {
      GraphBitset eGR2 = random_perm_graph(eGR);
      std::string str2 = get_string_expression(eGR2);
      if (str1 != str2) {
        std::cerr << "Error with the random permutation\n";
        throw TerminalException{1};
      }
    }
  }
  catch (TerminalException const& e) {
    exit(e.eVal);
  }
}
