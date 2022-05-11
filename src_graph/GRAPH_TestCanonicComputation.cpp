// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
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
    auto get_string_expression = [&](Tgr const &eGR) -> std::string {
      std::vector<Tidx> V;
      if (opt == 1) {
        std::cerr << "Running Traces canonicalization ordering\n";
        V = TRACES_GetCanonicalOrdering<Tgr, Tidx>(eGR);
      }
      if (opt == 2) {
        std::cerr << "Running Traces canonicalization ordering_arr_test\n";
        V = TRACES_GetCanonicalOrdering_Arr_Test<Tgr, Tidx>(eGR);
      }
      if (opt == 3) {
        std::cerr << "Running Bliss canonicalization ordering\n";
        V = BLISS_GetCanonicalOrdering<Tgr, Tidx>(eGR);
      }
      size_t nof_vertices = V.size();
      std::vector<unsigned int> clR(nof_vertices);
      for (size_t i = 0; i < nof_vertices; i++)
        clR[V[i]] = (unsigned int)i;
      std::string strRet;
      for (size_t iVert = 0; iVert < nof_vertices; iVert++) {
        int iVertCan = clR[iVert];
        strRet += " (";
        if (eGR.GetHasVertexColor()) {
          size_t eColor = eGR.GetColor(iVertCan);
          strRet += std::to_string(eColor) + ":";
        }
        //
        for (size_t jVert = 0; jVert < nof_vertices; jVert++) {
          size_t jVertCan = clR[jVert];
          bool eVal_b = eGR.IsAdjacent(iVertCan, jVertCan);
          strRet += " " + std::to_string(eVal_b);
        }
        strRet += ")";
      }
      return strRet;
    };
    auto random_perm_graph = [&](Tgr const &eGR,
                                 std::vector<Tidx> const &ePerm) -> Tgr {
      size_t nbVert = eGR.GetNbVert();
      Tgr eGR2(nbVert);
      bool HasVertexColor = eGR.GetHasVertexColor();
      if (HasVertexColor) {
        eGR2.SetHasColor(HasVertexColor);
        for (size_t iVert = 0; iVert < nbVert; iVert++) {
          size_t eColor = eGR2.GetColor(iVert);
          size_t iVert2 = size_t(ePerm[iVert]);
          eGR2.SetColor(iVert2, eColor);
        }
      }
      for (size_t iVert = 0; iVert < nbVert; iVert++) {
        size_t iVert2 = size_t(ePerm[iVert]);
        for (auto &eAdj : eGR.Adjacency(iVert)) {
          size_t eAdj2 = size_t(ePerm[eAdj]);
          eGR2.AddAdjacent(iVert2, eAdj2);
        }
      }
      return eGR2;
    };
    std::string str1 = get_string_expression(eGR);
    size_t nbVert = eGR.GetNbVert();
    for (int i = 0; i < 5; i++) {
      std::vector<Tidx> ePerm = RandomPermutation<Tidx>(nbVert);
      Tgr eGR2 = random_perm_graph(eGR, ePerm);
      std::string str2 = get_string_expression(eGR2);
      if (str1 != str2) {
        std::cerr << "str1=" << str1 << "\n";
        std::cerr << "str2=" << str2 << "\n";
        std::cerr << "Error with the random permutation\n";
        throw TerminalException{1};
      }
    }
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    exit(e.eVal);
  }
}
