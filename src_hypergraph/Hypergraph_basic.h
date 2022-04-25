#ifndef SRC_HYPERGRAPH_HYPERGRAPH_BASIC_H_
#define SRC_HYPERGRAPH_HYPERGRAPH_BASIC_H_

#include <vector>

struct HypergraphBitset {
  HypergraphBitset(int _nbVert, int _nbHyperedge,
                   std::vector<std::vector<int>> const &ListEdge)
      : nbVert(_nbVert), nbHyperedge(_nbHyperedge) {
    ListEntry = Face(nbVert * nbHyperedge);
    for (int iHyper = 0; iHyper < nbHyperedge; iHyper++) {
      for (auto &eVert : ListEdge[iHyper])
        ListEntry(iHyper * nbVert + eVert) = 1;
    }
  }
  int GetNbVert() const { return nbVert; }
  int GetNbHyperedge() const { return nbHyperedge; }
  std::vector<int> GetHyperedge(int const &iHyper) const {
    std::vector<int> LEnt;
    for (int iVert = 0; iVert < nbVert; iVert++) {
      int val = ListEntry(iHyper * nbVert + iVert);
      if (val == 1)
        LEnt.push_back(iVert);
    }
    return LEnt;
  }
  std::vector<int> GetContainingHyperedge(int const &iVert) const {
    std::vector<int> LEnt;
    for (int iHyper = 0; iHyper < nbHyperedge; iHyper++) {
      int val = ListEntry(iHyper * nbVert + iVert);
      if (val == 1)
        LEnt.push_back(iHyper);
    }
    return LEnt;
  }
  bool IsVertexInEdge(int const &iVert, int const &iHyper) const {
    return ListEntry(iHyper * nbVert + iVert);
  }

private:
  int nbVert;
  int nbHyperedge;
  Face ListEntry;
}

template <typename T, typename Tgr>
MyMatrix<T> GetLaplacianMatrix(std::vector<T> const &ListEdgeWeight,
                               Tgr const &GR) {
  int nbVert = GR.GetNbVert();
  std::vector<T> ListVertexWeight(nbVert);
  for (int iVert = 0; iVert < nbVert; iVert++) {
    std::vector<int> ListIEdge = GR.GetContainingHyperedge(iVert);
    T sum = 0;
    for (auto &iHyper : ListIEdge)
      sum += ListEdge[iHyper];
    ListVertexWeight[iVert] = sum;
  }
  //
  // Now building the matrix
  //
  int nbHyper = GR.GetNbHyperedge();
  MyMatrix<T> RetMat = ZeroMatrix<T>(nbVert, nbVert);
  for (int iHyper = 0; iHyper < nbHyper; iHyper++) {
    std::vector<int> ListIVert = GR.GetHyperedge(iHyper);
    int nbEnt = ListIVert.size();
    T coeff = ListEdgeWeight[iHyper] / nbEnt;
    for (int iEnt = 0; iEnt < nbEnt; iEnt++)
      for (int jEnt = iEnt + 1; jEnt < nbEnt; jEnt++) {
        int iVert = ListIVert[iEnt];
        int jVert = ListIVert[jEnt];
        T eWei = ListVertexWeight[iVert];
        T fWei = ListVertexWeight[jVert];
        RetMat(iVert, iVert) += coeff / eWei;
        RetMat(jVert, jVert) += coeff / fWei;
        T insCoeff = coeff / sqrt(eWei * fWei);
        RetMat(iVert, jVert) -= insCoeff;
        RetMat(jVert, iVert) -= insCoeff;
      }
  }
  return RetMat;
}

#endif // SRC_HYPERGRAPH_HYPERGRAPH_BASIC_H_
